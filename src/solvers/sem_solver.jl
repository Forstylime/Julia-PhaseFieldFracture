# src/solvers/sem_solver.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
    solve_sem(setup::TensionSetup, mat::PhaseFieldMaterial; kw...)

基于 S-EM (Self-Equilibrated Minimization) 算法的自适应时间步求解器。
使用弧长约束 ‖d - d_prev‖_M ≤ ρ 控制步长，能追踪 snap-back 响应。

三层嵌套循环:
- 外层: 自适应伪时间步进 (基于 L2 范数约束)
- 中层: 交错交替最小化 (位移 Newton → 相场约束求解)
- 内层: 增广拉格朗日 (Hestenes-Powell) + Sherman-Morrison 求解

通过 while 循环控制外层伪时间步, 直到 t_current 达到 t_max, 并记录 steps。

返回: (t_values, displacements, reaction_forces, elastic_energies, surface_energies)
"""
function solve_sem(
    setup::TensionSetup, mat::PhaseFieldMaterial;
    # n_steps::Int = 500,              # 最大自适应步数
    ρ::Float64 = 0.8658,               # 弧长半径
    α_init::Float64 = 10,           # AL 罚参数初始值
    β::Float64 = 1.2,                # 罚参数放大因子
    tol_staggered::Float64 = 1e-6,   # 交错收敛容差
    tol_al::Float64 = 1e-8,          # AL 内循环收敛容差
    max_staggered_iter::Int = 200,    # 最大交错迭代次数
    max_al_iter::Int = 200,           # 最大 AL 迭代次数
    max_newton_iter::Int = 10,       # 最大 Newton (位移) 迭代次数
    max_newton_iter_d::Int = 10,      # 最大 Newton (相场 d) 迭代次数
    tol_newton_d::Float64 = 1e-8,    # d-Newton 收敛容差
    t_max::Float64 = 86.58,            # 最大伪时间
    output_freq::Int = 5,            # VTK 输出频率
)
    # ================================================================
    # 1. 提取网格与自由度
    # ================================================================
    grid = setup.grid
    dh_u = setup.dh_u
    dh_d = setup.dh_d
    ch_u = setup.ch_u
    ch_d = setup.ch_d

    # 构建节点→位移自由度映射 (同 staggered.jl, 用于提取反力)
    node_dofs_u = zeros(Int, 2, getnnodes(grid))
    for cell_id in 1:getncells(grid)
        cell = getcells(grid, cell_id)
        dofs = celldofs(dh_u, cell_id)
        for (local_node, node_id) in pairs(cell.nodes)
            node_dofs_u[1, node_id] = dofs[(local_node - 1) * 2 + 1]
            node_dofs_u[2, node_id] = dofs[(local_node - 1) * 2 + 2]
        end
    end

    coords_x = [node.x[1] for node in grid.nodes]
    right_x = maximum(coords_x)
    right_nodes = findall(x -> isapprox(x, right_x; atol = 1e-12, rtol = 1e-12), coords_x)
    right_x_dofs = [node_dofs_u[2, node_id] for node_id in right_nodes]

    ndofs_u = ndofs(dh_u)
    ndofs_d = ndofs(dh_d)

    # ================================================================
    # 2. 初始化全局未知量
    # ================================================================
    u_n = zeros(ndofs_u);   u_prev = zeros(ndofs_u)
    d_n = zeros(ndofs_d);   d_prev = zeros(ndofs_d)

    apply!(u_n, ch_u)
    apply!(d_n, ch_d)
    u_prev .= u_n
    d_prev .= d_n

    # ================================================================
    # 3. 准备积分法则和 CellValues
    # ================================================================
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())

    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    H_history = zeros(n_qpoints)

    # ================================================================
    # 4. 分配稀疏矩阵和向量
    # ================================================================
    K_u = allocate_matrix(dh_u)
    K_d_base = allocate_matrix(dh_d)
    F_d_base = zeros(ndofs_d)
    R_u = zeros(ndofs_u)

    # 质量矩阵 (常数, 仅组装一次)
    M_d = allocate_matrix(dh_d)
    assemble_mass_matrix_d!(M_d, dh_d, cv_d)

    # Sherman-Morrison 工作空间
    A_sm  = allocate_matrix(dh_d)
    v_sm  = zeros(ndofs_d)
    r_mod = zeros(ndofs_d)
    x1    = zeros(ndofs_d)
    x2    = zeros(ndofs_d)
    dΔ    = zeros(ndofs_d)
    diff_d  = zeros(ndofs_d)
    M_diff  = zeros(ndofs_d)
    d_save  = zeros(ndofs_d)

    # ================================================================
    # 5. 记录结果用的数组
    # ================================================================
    t_values         = Float64[0.0]
    displacements    = Float64[0.0]
    reaction_forces  = Float64[0.0]
    elastic_energies  = Float64[0.0]
    surface_energies  = Float64[0.0]

    mkpath("data/sims")

    t_start = time()
    total_newton_iters = 0
    t_current = 0.0
    step = 0

    # ================================================================
    # 6. 外层: 自适应伪时间步循环
    # ================================================================
    println("开始 S-EM 求解, 最大步数: $n_steps, ρ = $ρ, t_max = $t_max")
    while t_current < t_max # for step in 1:n_steps # 
        # --- 6a. 时间步 ---
        step += 1
        if step == 1
            Δt = ρ
        else
            diff_d .= d_n .- d_prev
            mul!(M_diff, M_d, diff_d)
            Δd_norm = sqrt(abs(dot(diff_d, M_diff)))
            Δt = ρ - Δd_norm
        end

        t_current += Δt

        # 归一化伪时间用于 BC 缩放 (BC 期望 [0,1] 范围)
        normalized_t = t_current / t_max
        update!(ch_u, normalized_t)
        u_n .= u_prev
        d_n .= d_prev
        apply!(u_n, ch_u)
        apply!(d_n, ch_d)

        println("=== 步 $step | t = $(round(t_current, digits=6)) | Δt = $(round(Δt, digits=6)) ===")

        # ================================================================
        # 6b. 中层: 交错交替最小化循环
        # ================================================================
        staggered_converged = false
        for iter in 1:max_staggered_iter
            d_save .= d_n

            # ========================================
            # 步骤 A: 固定 d, Newton 求解位移 u
            # ========================================
            assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
            apply_zero!(K_u, R_u, ch_u)
            u_residual_norm = norm(R_u)
            newton_iter = 0

            while u_residual_norm > tol_staggered && newton_iter < max_newton_iter
                Δu = K_u \ (-R_u)
                apply_zero!(Δu, ch_u)
                u_n .+= Δu
                assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
                apply_zero!(K_u, R_u, ch_u)
                u_residual_norm = norm(R_u)
                newton_iter += 1
            end
            total_newton_iters += newton_iter

            # ========================================
            # 步骤 B: 更新不可逆历史变量 H
            # ========================================
            update_history!(H_history, dh_u, u_n, mat, cv_u)

            # ========================================
            # 步骤 C: 组装基础相场系统 (固定 H)
            # ========================================
            assemble_d!(K_d_base, F_d_base, dh_d, H_history, mat, cv_d)

            # ========================================
            # 步骤 D: AL 内循环求解带约束的 d
            # ========================================
            λ = 0.0        # 重置 Lagrange 乘子
            α = α_init     # 重置罚参数

            for hp_iter in 1:max_al_iter
                # ========================================
                # 内层: Newton 最小化 (λ, α 固定)
                # ========================================
                for newton_iter_d in 1:max_newton_iter_d
                    # 计算当前 d_k 下的约束量
                    diff_d .= d_n .- d_prev
                    mul!(M_diff, M_d, diff_d)
                    g = dot(diff_d, M_diff) - ρ^2

                    # 有效乘子 (使用外层的 λ, α, 不更新它们)
                    λ_eff = max(0.0, λ + α * g)

                    # 构造并求解 Newton 增量系统
                    if λ_eff > 0.0
                        # 约束 active: 使用 Sherman-Morrison 公式

                        # A = K_d_base + 2*λ_eff * M_d
                        @. A_sm.nzval = K_d_base.nzval + 2.0 * λ_eff * M_d.nzval

                        # ∇g = 2 * M_d * (d_n - d_prev) = 2 * M_diff
                        ∇g = 2.0 .* M_diff

                        # v = √α * ∇g
                        v_sm .= sqrt(α) .* ∇g

                        # r_mod = K_d_base*d_n - F_d_base + λ_eff*∇g
                        mul!(r_mod, K_d_base, d_n)
                        r_mod .-= F_d_base
                        r_mod .+= λ_eff .* ∇g

                        # 施加 Dirichlet BC
                        apply_zero!(A_sm, r_mod, ch_d)
                        apply_zero!(v_sm, ch_d)

                        # Sherman-Morrison: 一次分解求解两个右端项
                        b = -r_mod
                        sol = A_sm \ hcat(b, v_sm)
                        x1 .= sol[:, 1]
                        x2 .= sol[:, 2]

                        vTx1 = dot(v_sm, x1)
                        vTx2 = dot(v_sm, x2)

                        dΔ .= x1 .- (vTx1 / (1.0 + vTx2)) .* x2
                    else
                        # 约束 inactive: 标准线性求解
                        @. A_sm.nzval = K_d_base.nzval

                        mul!(r_mod, K_d_base, d_n)
                        r_mod .-= F_d_base

                        apply_zero!(A_sm, r_mod, ch_d)

                        dΔ .= A_sm \ (-r_mod)
                    end

                    # --- 更新 d ---
                    d_n .+= dΔ
                    d_n .= clamp.(max.(d_n, d_save), 0.0, 1.0)
                    apply!(d_n, ch_d)

                    # Newton 收敛检查
                    if norm(dΔ) < tol_newton_d
                        break
                    end
                end # Newton 内循环

                # ========================================
                # KKT 收敛检查 (Newton 收敛后)
                # ========================================
                diff_d .= d_n .- d_prev
                mul!(M_diff, M_d, diff_d)
                g = dot(diff_d, M_diff) - ρ^2
                viol = max(g, 0.0)
                comp = abs(λ * g)
                if viol < tol_al && comp < tol_al
                    break
                end

                # ========================================
                # Hestenes-Powell 乘子更新
                # ========================================
                λ = max(0.0, λ + α * g)
                α *= β
            end # HP 外循环

            # ========================================
            # 步骤 E: 检查交错收敛
            # ========================================
            d_error = norm(d_n - d_save)
            if mod(iter, 5) == 0
                @info " - 交错迭代 $iter, Δd_error = $(round(d_error, sigdigits=4))"
            end

            if d_error < tol_staggered
                staggered_converged = true
                break
            end
            if iter == max_staggered_iter
                @warn " - 步 $step 在最大交错迭代次数内未收敛！"
            end
        end # 中层交错循环

        # ================================================================
        # 6c. 保存状态
        # ================================================================
        u_prev .= u_n
        d_prev .= d_n

        # ================================================================
        # 6d. 后处理: 能量计算
        # ================================================================
        push!(elastic_energies, elastic_energy(dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d))
        push!(surface_energies, surface_energy(dh_d, d_n, mat, cv_d))

        # ================================================================
        # 6e. 反力计算
        # ================================================================
        assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
        f_total = sum(R_u[dof] for dof in right_x_dofs)

        push!(t_values, t_current)
        push!(displacements, normalized_t * setup.final_displacement)
        push!(reaction_forces, f_total)

        # ================================================================
        # 6f. VTK 输出
        # ================================================================
        if step % output_freq == 0
            VTKGridFile("data/sims/fracture_sem_step_$step", dh_u) do vtk
                write_solution(vtk, dh_u, u_n)
                write_solution(vtk, dh_d, d_n)
            end
        end

    end # 外层循环

    println("S-EM 仿真结束！VTK 文件保存在 data/sims/ 目录下。")
    println("总载荷步数: $step")
    println("总Newton迭代次数: $total_newton_iters")
    println("计算耗时: $(round(time() - t_start, digits=2)) 秒")

    return t_values, displacements, reaction_forces, elastic_energies, surface_energies
end
