# src/solvers/gamma.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
Gamma 弧长法求解器 (Monolithic)
具有双模态智能切换（弹性/断裂）与基于目标迭代次数的自适应步长控制。
"""
function solve_gamma(
    setup::MonolithicTensionSetup, mat::PhaseFieldMaterial;
    ρ_init::Float64 = 0.02,         # 初始断裂能步长
    Δλ_base_init::Float64 = 0.05,   # 初始位移步长（弹性阶段/尾盘阶段）
    max_steps::Int = 300,
    tol::Float64 = 1e-6,
    max_newton::Int = 15,
    output_freq::Int = 10,
    λ_max::Float64 = 1.0
)
    grid = setup.grid
    dh = setup.dh
    ch_ref = setup.ch_ref
    ch_zero = setup.ch_zero
    n_dofs = ndofs(dh)

    # ================================================================
    # 获取索引掩码
    # ================================================================
    idx_u = Int[]; idx_d = Int[]
    u_range = dof_range(dh, :u); d_range = dof_range(dh, :d)
    for cell in CellIterator(dh)
        cdofs = celldofs(cell)
        append!(idx_u, cdofs[u_range]); append!(idx_d, cdofs[d_range])
    end
    idx_u = unique(idx_u); idx_d = unique(idx_d)

    coords_x = [node.x[1] for node in grid.nodes]; right_x = maximum(coords_x)
    right_dofs = Int[]
    for cell in CellIterator(dh)
        for (i, node_id) in enumerate(cell.nodes)
            if isapprox(grid.nodes[node_id].x[1], right_x; atol = 1e-8)
                push!(right_dofs, celldofs(cell)[u_range[(i-1)*2 + 2]])
            end
        end
    end
    right_dofs = unique(right_dofs)

    # 状态初始化 
    a_n = zeros(n_dofs); a_cur = zeros(n_dofs)
    λ_n = 0.0; λ_cur = 0.0
    G_prev = 0.0  

    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    H_old = zeros(n_qpoints)

    K_mono = allocate_matrix(dh); r_mono = zeros(n_dofs)
    
    displacements = Float64[0.0]; reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]; surface_energies = Float64[0.0]

    mkpath("data/sims/gamma")
    
    # 【全局控制变量初始化】
    ρ = ρ_init
    Δλ_base = Δλ_base_init
    n_step = 1
    t_start = time()

    println("开始 Gamma 弧长法，ρ_init = $ρ_init, Δλ_init = $Δλ_base_init")

    while λ_n <= λ_max && n_step < max_steps
        println("=== 载荷步 $n_step | 目标 λ ≈ $(round(λ_n + Δλ_base, digits=4)) ===")

        # =============================================
        # PHASE 1: PREDICTOR (预测步)
        # =============================================
        assemble_monolithic!(K_mono, r_mono, dh, a_n, H_old, mat, cv_u, cv_d)
        K_f_pred = copy(K_mono); f_f_pred = zeros(n_dofs)
        apply!(K_f_pred, f_f_pred, ch_ref)
        δa_λ_pred = K_f_pred \ f_f_pred

        _, K_λa_n, _ = evaluate_gamma_constraint(dh, a_n, G_prev, mat, cv_d, 0.0)
        sens_G_λ = dot(K_λa_n, δa_λ_pred)
        
        # 智能判定：当系统几乎没有能量，或者断裂能对位移脱敏时，进入位移控制
        is_gamma_active = (G_prev > 1e-3) && (abs(sens_G_λ) > 1e-2)

        if !is_gamma_active
            δλ_pred = Δλ_base
            println("  -> [模式: 弹性/尾盘] 预测 Δλ = $(round(δλ_pred, digits=4))")
        else
            δλ_pred = ρ / sens_G_λ
            # 弧长法预测器限幅：防止刚开裂时一脚油门冲太猛
            δλ_pred = clamp(δλ_pred, -0.1, 0.1) 
            println("  -> [模式: Γ 弧长] 预测 Δλ = $(round(δλ_pred, digits=4)), 当前 ρ = $(round(ρ, digits=4))")
        end
        
        δa_pred = δλ_pred .* δa_λ_pred 
        λ_cur = λ_n + δλ_pred
        a_cur .= a_n .+ δa_pred

        # =============================================
        # PHASE 2: CORRECTOR (Full Newton 校正步)
        # =============================================
        converged = false
        iters_newton = 0
        for iter in 1:max_newton
            iters_newton = iter
            assemble_monolithic!(K_mono, r_mono, dh, a_cur, H_old, mat, cv_u, cv_d)

            r_check = -copy(r_mono); apply_zero!(r_check, ch_zero)
            if norm(r_check) <= tol; converged = true; break; end

            K_λ = copy(K_mono); f_λ = zeros(n_dofs)
            apply!(K_λ, f_λ, ch_ref); δa_λ = K_λ \ f_λ

            K_r = copy(K_mono); r = -copy(r_mono)
            apply!(K_r, r, ch_zero); δa_r = K_r \ r

            # 在校正步中严格执行对应的约束，无 clamp 干预，确保二次收敛
            if !is_gamma_active
                δλ = 0.0
            else
                f_Γ, K_λa_iter, _ = evaluate_gamma_constraint(dh, a_cur, G_prev, mat, cv_d, ρ)
                denominator = dot(K_λa_iter, δa_λ)
                
                if abs(denominator) < 1e-14
                    δλ = 0.0
                else
                    numerator = f_Γ + dot(K_λa_iter, δa_r)
                    δλ = - numerator / denominator 
                end
            end
           
            δa_total = δa_r .+ δλ .* δa_λ
            a_cur .+= δa_total
            λ_cur += δλ
        end 

        # =============================================
        # PHASE 3: 收敛后处理 & 动态自适应步长
        # =============================================
        if converged
            # 收敛后结算这一步真实的断裂能
            _, _, G_curr = evaluate_gamma_constraint(dh, a_cur, 0.0, mat, cv_d, 0.0)
            G_prev = G_curr
            
            # 目标迭代缩放法则 (工业界标配)
            target_iters = 4.0
            factor = sqrt(target_iters / max(iters_newton, 1))
            factor = clamp(factor, 0.5, 2.0) # 限制单次突变不超过2倍
            
            if is_gamma_active
                ρ = ρ * factor
                ρ = min(ρ, 5.0) # 物理上限：L-shape 断裂能总共80，每步走5.0非常快
            else
                Δλ_base = Δλ_base * factor
                Δλ_base = clamp(Δλ_base, 0.01, 0.1) # 弹性区最大一步走 0.1
            end
            
            println("  -> 成功收敛 (用时 $iters_newton 步). 最终 λ = $(round(λ_cur, digits=4))")
            
            update_history_mono!(H_old, dh, a_cur, mat, cv_u)
            a_n .= a_cur
            λ_n = λ_cur

            assemble_monolithic!(K_mono, r_mono, dh, a_n, H_old, mat, cv_u, cv_d)
            f_reac = sum(r_mono[dof] for dof in right_dofs)

            push!(displacements, λ_n * setup.final_displacement)
            push!(reaction_forces, f_reac)
            push!(elastic_energies, elastic_energy_monolithic(dh, a_n, mat, cv_u, cv_d))
            push!(surface_energies, G_prev)

            if n_step % output_freq == 0
                VTKGridFile("data/sims/gamma/gamma_step_$n_step", setup.grid) do vtk
                    write_solution(vtk, dh, a_n) 
                end
            end
            n_step += 1
        else
            @warn "未在 $max_newton 次迭代内收敛。回退状态并减小步长..."
            # 针对不同模式安全缩减步长
            if is_gamma_active
                ρ *= 0.5
                if ρ < 1e-6; error("弧长 ρ 过小，仿真发散中止。"); end
            else
                Δλ_base *= 0.5
                if Δλ_base < 1e-5; error("弹性位移步长过小，仿真发散中止。"); end
            end
            continue # 跳过状态更新，重新在 a_n 的基础上执行 Predictor
        end
    end

    println("Gamma 仿真结束！计算耗时: $(round(time() - t_start, digits=2)) 秒, 总步数: $(n_step - 1)")
    return displacements, reaction_forces, elastic_energies, surface_energies
end