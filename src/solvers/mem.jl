# src/solvers/mem.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
    solve_mem(setup::MonolithicTensionSetup_MEM, mat::PhaseFieldMaterial; kw...)

基于 M-EM (Monolithic Efendiev & Mielke) 算法的自适应时间步求解器。
"""
function solve_mem(
    setup::MonolithicTensionSetup_MEM, mat::PhaseFieldMaterial;
    max_steps::Int = 500,
    ρ_init::Float64 = 0.8658,             # 相场演化上限 (弧长)
    tol_newton::Float64 = 1e-4,      # Newton 迭代收敛容差
    tol_kkt::Float64 = 1e-6,         # 状态切换判断容差
    max_newton::Int = 30,            # 最大 Newton 迭代次数
    t_max::Float64 = 86.58,          # 最大伪时间 (对应最终位移)
    output_freq::Int = 10,
)
    # ================================================================
    # 1. 基础网格与自由度提取 
    # ================================================================
    grid = setup.grid
    dh = setup.dh
    
    ch_a = setup.ch_a 
    ch_zero = setup.ch_zero # 用于 Newton 更新步的齐次 BC
    n_dofs = ndofs(dh)

    # 提取 u 和 d 在交错全局向量中的索引
    idx_u = Int[]; idx_d = Int[]
    u_range = dof_range(dh, :u); d_range = dof_range(dh, :d)
    for cell in CellIterator(dh)
        cdofs = celldofs(cell)
        append!(idx_u, cdofs[u_range])
        append!(idx_d, cdofs[d_range])
    end
    sort!(unique!(idx_u)); sort!(unique!(idx_d))

    # 提取右端点位移自由度（用于计算反力）
    coords_x = [node.x[1] for node in grid.nodes]
    right_x = maximum(coords_x)
    right_dofs = Int[]
    for cell in CellIterator(dh)
        for (i, node_id) in enumerate(cell.nodes)
            if isapprox(grid.nodes[node_id].x[1], right_x; atol = 1e-8)
                push!(right_dofs, celldofs(cell)[u_range[(i-1)*2 + 2]])
            end
        end
    end
    right_dofs = unique(right_dofs)

    # ================================================================
    # 2. 积分与历史变量初始化
    # ================================================================
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())

    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    H_old = zeros(n_qpoints)

    # 提取单场的 dh_d 组装纯相场的质量矩阵 M_d
    dh_d = setup.dh_d 
    M_d = allocate_matrix(dh_d)
    assemble_mass_matrix_d!(M_d, dh_d, cv_d)

    # 全局矩阵与向量
    K_mono = allocate_matrix(dh)
    r_mono = zeros(n_dofs)

    # ================================================================
    # 3. 状态回溯数组初始化 (M-EM 核心)
    # ================================================================
    a_prev_prev = zeros(n_dofs)
    a_prev      = zeros(n_dofs)
    a_cur       = zeros(n_dofs)

    # 零时刻边界条件初始化
    update!(ch_a, 0.0)
    apply!(a_prev, ch_a)
    a_prev_prev .= a_prev

    displacements = Float64[0.0]
    reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]
    surface_energies = Float64[0.0]

    mkpath("data/sims/mem")
    t_start = time()
    t_prev = 0.0; t_cur = 0.0
    n_step = 1
    total_newton_iters = 0
    ρ = ρ_init

    println("开始 M-EM 求解, 初始 ρ_0 = $ρ_init, t_max = $t_max")

    # ================================================================
    # 4. 外层时间自适应循环 (Algorithm 2)
    # ================================================================
    while t_cur < t_max && n_step <= max_steps
        println("\n=== 载荷步 $n_step: t_cur = $(round(t_cur, digits=4)) ===")
        
        step_converged = false

        # 内层尝试循环：若当前步不收敛，减小 ρ 后重新计算 t_cur 并重试
        while !step_converged
            # --- 4.1 物理伪时间更新 (Eq 38c) ---
            if n_step == 1
                t_cur = ρ
            else
                Δd_old = a_prev[idx_d] .- a_prev_prev[idx_d]
                norm_d_old = sqrt(dot(Δd_old, M_d * Δd_old))
                t_cur = min(t_prev + ρ - norm_d_old, t_max)
            end

            # --- 4.2 更新并强加 Dirichlet 边界条件 ---
            normalized_t = t_cur / t_max
            update!(ch_a, normalized_t) 

            # --- 4.3 构造当前步的预测器 Predictor (Eq 48, 49) ---
            if n_step == 1
                a_cur .= a_prev
            else
                a_cur .= 2.0 .* a_prev .- a_prev_prev
            end
            apply!(a_cur, ch_a) # 将当前时间的边界条件强加给预测器

            # ================================================================
            # 5. 内层：有效集与牛顿求解循环
            # ================================================================
            active_converged = false
            λ_rho = 0.0
            
            # 初始激活状态猜测
            Δd_pred = a_cur[idx_d] .- a_prev[idx_d]
            g_val = dot(Δd_pred, M_d * Δd_pred) - ρ^2
            is_active = g_val >= 0.0 

            newton_failed = false
            active_attempts = 0 # 用于防止 Active-Inactive 状态切换振荡

            while !active_converged
                active_attempts += 1
                if active_attempts > 10
                    @warn "    [警告] 有效集状态切换频繁，可能发生非凸振荡，强制触发回溯..."
                    newton_failed = true
                    break
                end

                if is_active
                    println("  -> [约束 Active]: 剧烈扩展阶段，开启弧长控制...")
                    converged, a_cur, λ_rho, iters = solve_newton_active!(
                        a_cur, a_prev, λ_rho, K_mono, r_mono, dh, ch_zero, ch_a, H_old, mat, cv_u, cv_d, 
                        M_d, ρ, idx_u, idx_d, max_newton, tol_newton
                    )
                    total_newton_iters += iters

                    if !converged
                        newton_failed = true
                        break 
                    end

                    # KKT 检查 (仅检查拉格朗日乘子正负，去掉无意义的 g_val 校验)
                    if λ_rho < -tol_kkt
                        println("     [KKT 检查] λ = $(round(λ_rho, digits=6)) < -$tol_kkt. 退回为 Inactive。")
                        is_active = false
                    else
                        active_converged = true
                    end
                else
                    println("  -> [约束 Inactive]: 平稳扩展阶段，标准 Newton 求解...")
                    λ_rho = 0.0
                    converged, a_cur, iters = solve_newton_inactive!(
                        a_cur, K_mono, r_mono, dh, ch_zero, H_old, mat, cv_u, cv_d, max_newton, tol_newton
                    )
                    total_newton_iters += iters

                    if !converged
                        # 【重要逻辑修正】：若标准求解发散，说明发生非稳态跳跃，不回溯，直接切入 Active 求解
                        println("     [求解提示] Inactive 求解发散，可能发生非稳态跳跃，尝试切换至 Active 求解...")
                        is_active = true
                        continue
                    end

                    # KKT 检查
                    Δd_res = a_cur[idx_d] .- a_prev[idx_d]
                    g_res = dot(Δd_res, M_d * Δd_res) - ρ^2
                    if g_res >= tol_kkt
                        println("     [KKT 检查] g(d) = $(round(g_res, digits=4)) >= $tol_kkt. 演化超限，切入 Active。")
                        is_active = true
                    else
                        active_converged = true
                    end
                end
            end # 有效集循环结束

            # 处理真正的不收敛情况 (例如 Active 求解因收敛半径问题失败)
            if newton_failed
                @warn "载荷步 $n_step Newton 求解失败，自动减小 ρ 并完全回溯重试该步"
                ρ *= 0.5
                if ρ < 1e-6
                    error("弧长 ρ 过小 ($(ρ))，计算终止。")
                end
            else
                ρ = ρ_init # 适当增加 ρ 以加速后续步长
                step_converged = true # 成功求解，跳出 retry 循环
            end

        end # retry 循环结束

        # ================================================================
        # 6. 后处理与状态更新
        # ================================================================
        update_history_mono!(H_old, dh, a_cur, mat, cv_u)

        a_prev_prev .= a_prev
        a_prev .= a_cur
        t_prev = t_cur

        # 重新装配以计算平衡状态下的反力与能量
        assemble_monolithic!(K_mono, r_mono, dh, a_cur, H_old, mat, cv_u, cv_d)
        f_reac = sum(r_mono[dof] for dof in right_dofs)

        normalized_t = t_cur / t_max
        push!(displacements, normalized_t * setup.final_displacement)
        push!(reaction_forces, f_reac)
        push!(elastic_energies, elastic_energy_monolithic(dh, a_cur, mat, cv_u, cv_d))
        push!(surface_energies, surface_energy_monolithic(dh, a_cur, mat, cv_d))

        if n_step % output_freq == 0
            VTKGridFile("data/sims/mem/mem_step_$n_step", setup.grid) do vtk
                write_solution(vtk, dh, a_cur)
            end
        end
        n_step += 1
    end

    println("M-EM 仿真结束！计算耗时: $(round(time() - t_start, digits=2)) 秒, 总 Newton 迭代: $total_newton_iters")
    return displacements, reaction_forces, elastic_energies, surface_energies
end