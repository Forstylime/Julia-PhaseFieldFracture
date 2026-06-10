# src/solvers/L2.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
L2 弧长法求解器 (Monolithic)
直接采用纯线性代数操作(M_mat)评估约束，无需在牛顿迭代中重算约束雅可比。
"""
function solve_l2(
    setup::MonolithicTensionSetup, mat::PhaseFieldMaterial;
    # 初始步长
    ρ_init::Float64 = 0.05,         
    Δλ_base_init::Float64 = 0.05,   
    
    # ======= 新增：精细度控制阀门 =======
    ρ_max::Float64 = 5.0,          # (核心) 断裂期间的最大弧长限制。越小，裂纹扩展阶段的步数越多、曲线越平滑
    Δλ_max::Float64 = 0.2,         # 弹性/未开裂期间的最大位移步长。越小，弹性直线上点越多
    target_iters::Float64 = 5.0,   # 自适应激进程度。设为 3.0 会更趋向于保守切分
    # ====================================

    max_steps::Int = 1000,         # 注意把最大步数放宽，防止高精模式下跑不完

    # 以下是常规参数，保持不变
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

    # 获取右侧边界自由度等逻辑与 Gamma 法完全一致
    coords_x = [node.x[1] for node in grid.nodes]; right_x = maximum(coords_x)
    right_dofs = Int[]
    for cell in CellIterator(dh)
        for (i, node_id) in enumerate(cell.nodes)
            if isapprox(grid.nodes[node_id].x[1], right_x; atol = 1e-8)
                u_range = dof_range(dh, :u)
                push!(right_dofs, celldofs(cell)[u_range[(i-1)*2 + 2]])
            end
        end
    end
    right_dofs = unique(right_dofs)

    # 状态初始化 
    a_n = zeros(n_dofs); a_cur = zeros(n_dofs)
    λ_n = 0.0; λ_cur = 0.0

    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    H_old = zeros(n_qpoints)

    K_mono = allocate_matrix(dh); r_mono = zeros(n_dofs)
    
    # 【核心修改 1：在循环外一次性装配全局 H 矩阵】
    println("正在预计算恒定 L2 质量矩阵...")
    M_mat = assemble_L2_matrix(dh, cv_d)
    
    displacements = Float64[0.0]; reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]; surface_energies = Float64[0.0]

    mkpath("data/sims/l2")
    
    ρ = ρ_init
    Δλ_base = Δλ_base_init
    n_step = 1
    t_start = time()
    has_fractured_before = false

    println("开始 L2 弧长法，ρ_init = $ρ_init, Δλ_init = $Δλ_base_init")

    while λ_n <= λ_max && n_step < max_steps
        println("=== 载荷步 $n_step | 目标 λ ≈ $(round(λ_n + Δλ_base, digits=4)) ===")

        # =============================================
        # PHASE 1: PREDICTOR (预测步)
        # =============================================
        assemble_monolithic!(K_mono, r_mono, dh, a_n, H_old, mat, cv_u, cv_d)
        K_f_pred = copy(K_mono); f_f_pred = zeros(n_dofs)
        apply!(K_f_pred, f_f_pred, ch_ref)
        δa_λ_pred = K_f_pred \ f_f_pred

        # 评估当前预测的相场变化趋势 (D_pred = Δd^T * M * Δd)
        D_pred = dot(δa_λ_pred, M_mat * δa_λ_pred)
        
        # 智能判定：(D_pred > 1e-10) 且有实质性载荷
        is_l2_active = (D_pred > 1e-10) && (λ_n > 0.1)

        if !is_l2_active
            δλ_pred = Δλ_base
            println("  -> [模式: 弹性] 预测 Δλ = $(round(δλ_pred, digits=4))")
        else
            # 【关键优化】智能无缝切换：第一次进入断裂时，动态继承弹性的速度
            if !has_fractured_before
                # 用弹性的最后一次预测步长，反推一个等效的初始 ρ
                ρ = abs(Δλ_base * sqrt(D_pred))
                has_fractured_before = true
                println("  -> [状态转换] 首次起裂！动态匹配初始 ρ = $(round(ρ, digits=4))")
            end

            sum_d_pred = sum(δa_λ_pred[idx_d])
            sign_λ = sum_d_pred >= 0.0 ? 1.0 : -1.0 
            
            δλ_pred = sign_λ * (ρ / sqrt(D_pred))
            
            # 【解除预测步限幅】原先是 clamp(-0.1, 0.1)，或许太保守了！
            # 可以考虑允许弧长法在预测时给出更大的步长
            δλ_pred = clamp(δλ_pred, -0.1, 0.1) 
            println("  -> [模式: L2 弧长] 预测 Δλ = $(round(δλ_pred, digits=4)) (Snap-back 标志: $sign_λ)")
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

            if !is_l2_active
                δλ = 0.0
            else
                # 【修复点 2：采用非平方 Norm 形式，彻底解决梯度归零/奇异性】
                Δa_cur = a_cur .- a_n
                D_cur = dot(Δa_cur, M_mat * Δa_cur)
                
                if D_cur < 1e-15
                    # 如果增量极小，强制降级为纯位移校正，避免 0/0 错位
                    δλ = 0.0
                else
                    norm_L2 = sqrt(D_cur)
                    f_L2_val = norm_L2 - ρ
                    
                    # 单位化梯度： K_λa_iter = (M * Δa) / norm_L2
                    K_λa_iter = (M_mat * Δa_cur) ./ norm_L2
                    
                    denominator = dot(K_λa_iter, δa_λ)
                    
                    # 【修复点 3：分母截断保护机制】
                    if abs(denominator) < 1e-12
                        @warn "  Newton iter $iter: 切线分母趋近零，临时冻结载荷步"
                        δλ = 0.0
                    else
                        numerator = f_L2_val + dot(K_λa_iter, δa_r)
                        δλ = - numerator / denominator 

                        δλ = clamp(δλ, -0.1, 0.1)
                    end
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
            # 使用传入的 target_iters 来计算缩放因子
            factor = sqrt(target_iters / max(iters_newton, 1))
            factor = clamp(factor, 0.25, 2.5) 
            
            if is_l2_active
                # 【核心控制 1：断裂期精细度】
                # 利用 ρ_max 牢牢压住弧长上限。如果 ρ_max 设得很小（比如 0.1），
                # 哪怕牛顿法 1 步就收敛，系统也只能乖乖地一次走 0.1 的弧长
                ρ = clamp(ρ * factor, 1e-4, ρ_max) 
            else
                # 【核心控制 2：弹性期精细度】
                # 利用 Δλ_max 压住载荷加载速度
                Δλ_base = clamp(Δλ_base * factor, 0.01, Δλ_max)
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
            push!(surface_energies, surface_energy_monolithic(dh, a_n, mat, cv_d)) # 记录断裂能

            if n_step % output_freq == 0
                VTKGridFile("data/sims/l2/l2_step_$n_step", setup.grid) do vtk
                    write_solution(vtk, dh, a_n) 
                end
            end
            n_step += 1
        else
            @warn "未在 $max_newton 次迭代内收敛。回退状态并减小步长..."
            if is_l2_active
                ρ *= 0.5
                if ρ < 1e-6; error("弧长 ρ 过小，仿真发散中止。"); end
            else
                Δλ_base *= 0.5
                if Δλ_base < 1e-6; error("位移步长过小，仿真发散中止。"); end
            end
            continue
        end
    end

    println("L2 仿真结束！计算耗时: $(round(time() - t_start, digits=2)) 秒, 总步数: $(n_step - 1)")
    return displacements, reaction_forces, elastic_energies, surface_energies
end