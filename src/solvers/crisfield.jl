# src/solvers/crisfield.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
Crisfield 弧长法求解器
"""
function solve_crisfield(
    setup::MonolithicTensionSetup, mat::PhaseFieldMaterial;
    ρ_init::Float64 = 0.01,
    max_steps::Int = 500,
    tol::Float64 = 1e-6,
    max_newton::Int = 15,
    output_freq::Int = 10,
    λ_max::Float64 = 1.0,
    enforce_irreversibility::Bool = true
)
    dir = setup.dir
    grid = setup.grid
    dh = setup.dh
    ch_ref = setup.ch_ref
    ch_zero = setup.ch_zero
    n_dofs = ndofs(dh)

    # ================================================================
    # 获取 u 和 d 在全局向量中的索引掩码 
    # ================================================================
    idx_u = Int[]
    idx_d = Int[]
    u_range = dof_range(dh, :u)
    d_range = dof_range(dh, :d)
    for cell in CellIterator(dh)
        cdofs = celldofs(cell)
        append!(idx_u, cdofs[u_range])
        append!(idx_d, cdofs[d_range])
    end
    idx_u = unique(idx_u); idx_d = unique(idx_d)

    coords_x = [node.x[1] for node in grid.nodes]
    right_x = maximum(coords_x)
    right_dofs = Int[]
    for cell in CellIterator(dh)
        for (i, node_id) in enumerate(cell.nodes)
            if isapprox(grid.nodes[node_id].x[1], right_x; atol = 1e-8)
                push!(right_dofs, celldofs(cell)[u_range[(i-1)*2 + dir]])
            end
        end
    end
    right_dofs = unique(right_dofs)

    # 初始化状态 
    a_n = zeros(n_dofs)
    a_cur = zeros(n_dofs)
    λ_n = 0.0; λ_cur = 0.0

    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    driving_force = zeros(n_qpoints)

    K_mono = allocate_matrix(dh)
    r_mono = zeros(n_dofs)
    
    displacements = Float64[0.0]
    reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]
    surface_energies = Float64[0.0]

    mkpath("data/sims/crisfield2")
    ρ = ρ_init
    Δa_n = zeros(n_dofs) # 上一个收敛步的总增量，用于判断前进方向
    n_step = 1
    t_start = time()

    println("开始 Crisfield 弧长法 (Monolithic)，ρ_init = $ρ, λ_max = $λ_max")

    while λ_n <= λ_max && n_step < max_steps
        println("=== 载荷步 $n_step | λ = $(round(λ_n, digits=4)) | ρ = $(round(ρ, digits=4)) ===")

        # =============================================
        # PHASE 1: PREDICTOR (预测步)
        # =============================================
        assemble_monolithic!(K_mono, r_mono, dh, a_n, driving_force, mat, cv_u, cv_d)

        K_f_pred = copy(K_mono)
        f_f_pred = zeros(n_dofs)
        apply!(K_f_pred, f_f_pred, ch_ref)
        
        δa_λ_pred = K_f_pred \ f_f_pred

        # 弧长缩放必须只针对位移自由度！
        norm_δu_λ_pred = norm(δa_λ_pred[idx_u])
        abs_δλ_pred = ρ / norm_δu_λ_pred

        # 方向判断点积也只针对位移自由度
        if n_step == 1
            sign_δλ = 1.0
        else
            sign_δλ = dot(Δa_n[idx_u], δa_λ_pred[idx_u]) + dot(Δa_n[idx_d], δa_λ_pred[idx_d]) >= 0.0 ? 1.0 : -1.0
        end
        
        δλ_pred = sign_δλ * abs_δλ_pred
        δa_pred = δλ_pred .* δa_λ_pred # 乘全局向量，带动 d 场猜测

        # 初始猜测状态 
        λ_cur = λ_n + δλ_pred
        a_cur .= a_n .+ δa_pred
        Δa_iter = copy(δa_pred) # 当前增量步的累计增量

        # =============================================
        # PHASE 2: CORRECTOR (Full Newton)
        # =============================================
        converged = false
        iters_newton = 0
        for iter in 1:max_newton
            iters_newton = iter
            assemble_monolithic!(K_mono, r_mono, dh, a_cur, driving_force, mat, cv_u, cv_d)

            # --- 2.1. 收敛检查 ---
            r_check = -copy(r_mono)
            apply_zero!(r_check, ch_zero)
            res_norm = norm(r_check)

            if res_norm <= tol
                converged = true
                break
            end

            # --- 2.2. 计算载荷方向增量 δa_λ ---
            K_λ = copy(K_mono)
            f_λ = zeros(n_dofs)
            apply!(K_λ, f_λ, ch_ref)
            δa_λ = K_λ \ f_λ

            # --- 2.3. 计算残差方向增量 δa_r ---
            K_r = copy(K_mono)
            r = -copy(r_mono)
            apply!(K_r, r, ch_zero)
            δa_r = K_r \ r

            # --- 2.4. 采用一致线性化公式 (26) 求解标量 δλ (只针对位移自由度 idx_u) ---
            Δu_iter = Δa_iter[idx_u]
            δu_r = δa_r[idx_u]
            δu_λ = δa_λ[idx_u]

            # 计算 1: 当前弧长约束的偏差标量值 f_bullet
            f_bullet = dot(Δu_iter, Δu_iter) - ρ^2

            # 计算 2: 导数与增量的点积
            # 因为约束 f = Δu_iter^2 - ρ^2 对 u 的梯度是 2 * Δu_iter
            # 所以 K_la * δa_r 即为 2 * dot(Δu_iter, δu_r)
            # K_la * δa_λ 即为 2 * dot(Δu_iter, δu_λ)
            # 此外，K_ll = 0.0 (因为约束中不包含载荷因子 λ)
            
            numerator = -(f_bullet + 2.0 * dot(Δu_iter, δu_r))
            denominator = 2.0 * dot(Δu_iter, δu_λ)

            # 3. 得到线性化步长 δλ
            δλ = numerator / denominator
           
            # --- 2.5. 更新状态 ---
            δa_total = δa_r .+ δλ .* δa_λ
            Δa_iter .+= δa_total
            a_cur .+= δa_total
            λ_cur += δλ
        end # Newton 迭代结束

        # =============================================
        # PHASE 3: 收敛后处理 → 步长控制
        # =============================================
        if converged
            Δa_n .= Δa_iter # 记录本步收敛的纯增量，供下一步判断方向
            
            # 动态放大弧长，加速收敛
            if iters_newton <= 3
                ρ = min(ρ * 1.2, ρ_init * 5.0) # 设一个最大弧长限制
            end
            if n_step == 150
                ρ_init = ρ_init / 10
            end
        else
            @warn "未在 $max_newton 次迭代内收敛。回退并减半弧长 ρ。"
            ρ *= 0.5
            continue # a_n 没被污染，直接重做 Predictor
        end

        # 更新历史变量, driving_force 代表论文力的 H
        compute_driving_force_mono!(driving_force, dh, a_cur, mat, cv_u, enforce_irreversibility)

        a_n .= a_cur
        λ_n = λ_cur

        # 为了计算严谨的反力和能量，最后用收敛态装配一次
        assemble_monolithic!(K_mono, r_mono, dh, a_n, driving_force, mat, cv_u, cv_d)
        f_reac = sum(r_mono[dof] for dof in right_dofs)

        # 记录数据 (这里 final_displacement 是参考基准)
        push!(displacements, λ_n * setup.final_displacement)
        push!(reaction_forces, f_reac)
        push!(elastic_energies, elastic_energy_monolithic(dh, a_n, mat, cv_u, cv_d))
        push!(surface_energies, surface_energy_monolithic(dh, a_n, mat, cv_d)) 

        if n_step % output_freq == 0
            VTKGridFile("data/sims/crisfield2/crisfield_step_$n_step", setup.grid) do vtk
                write_solution(vtk, dh, a_n) 
            end
        end
        n_step += 1
    end

    println("Crisfield 仿真结束！VTK 文件保存在 data/sims/crisfield2 目录下。")
    println("计算耗时: $(round(time() - t_start, digits=2)) 秒")
    return displacements, reaction_forces, elastic_energies, surface_energies
end