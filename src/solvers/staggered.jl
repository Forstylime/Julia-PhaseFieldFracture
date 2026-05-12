# src/solvers/staggered.jl

using Ferrite
using LinearAlgebra
using SparseArrays

function _nodal_dof_map(dh::DofHandler; components::Int = 1)
    grid = Ferrite.get_grid(dh)
    node_dofs = zeros(Int, components, getnnodes(grid))

    for cell_id in 1:getncells(grid)
        cell = getcells(grid, cell_id)
        dofs = celldofs(dh, cell_id)
        for (local_node, node_id) in pairs(cell.nodes)
            for component in 1:components
                node_dofs[component, node_id] = dofs[(local_node - 1) * components + component]
            end
        end
    end

    return node_dofs
end

function _crack_constraint_handler(dh_d::DofHandler, crack_nodes)
    ch_d = ConstraintHandler(dh_d)
    if !isempty(crack_nodes)
        add!(ch_d, Dirichlet(:d, collect(crack_nodes), (x, t) -> 1.0))
    end
    close!(ch_d)
    update!(ch_d, 0.0)
    return ch_d
end

function _top_component_dofs(dh_u::DofHandler, grid, component::Int)
    coords_y = [node.x[2] for node in grid.nodes]
    top_y = maximum(coords_y)
    top_nodes = findall(y -> isapprox(y, top_y; atol = 1e-12, rtol = 1e-12), coords_y)
    node_dofs = _nodal_dof_map(dh_u; components = 2)
    return [node_dofs[component, node_id] for node_id in top_nodes]
end

"""
    solve_staggered(setup::SquareTensionSetup, mat::PhaseFieldMaterial; 
                    n_steps = 100, max_u_disp = 0.01, tol = 1e-4)

执行标准交错求解法 (Staggered Scheme)。
在每个载荷步中，交错求解位移场和相场，直到残差收敛。
"""
function solve_staggered(
    setup::SquareTensionSetup, mat::PhaseFieldMaterial;
    n_steps = 100,          # 总载荷步数 (论文设为 100 步)
    max_u_disp = 0.01,      # 顶部的最大位移加载量 (mm)
    tol = 1e-4,             # Newton 迭代和交错循环的收敛容差
    max_iter = 20           # 内层交错循环的最大允许次数
)
    # --- 1. 提取网格与自由度 ---
    grid = setup.grid
    dh_u = setup.dh_u
    dh_d = setup.dh_d
    ch_u = create_displacement_constraints(dh_u, grid; top_displacement = max_u_disp)
    ch_d = _crack_constraint_handler(dh_d, setup.crack_nodes)
    top_y_dofs = _top_component_dofs(dh_u, grid, 2)
    
    ndofs_u = ndofs(dh_u)
    ndofs_d = ndofs(dh_d)
    
    # 初始化全局未知量
    u_n = zeros(ndofs_u); u_prev = zeros(ndofs_u)
    d_n = zeros(ndofs_d); d_prev = zeros(ndofs_d)
    
    # 初始化边界条件与预制裂纹。裂纹节点通过 Dirichlet 条件固定为 d = 1。
    apply!(u_n, ch_u)
    apply!(d_n, ch_d) 
    u_prev .= u_n
    d_prev .= d_n
    
    # --- 2. 准备积分法则和 CellValues ---
    # 四边形单元，2阶高斯积分 (2x2=4个积分点)，对于相场断裂精度足够
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    
    # 计算总积分点数量，初始化历史变量 H
    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    H_history = zeros(n_qpoints)
    
    # 分配全局稀疏矩阵和向量
    K_u = allocate_matrix(dh_u)
    K_d = allocate_matrix(dh_d)
    R_u = zeros(ndofs_u)
    F_d = zeros(ndofs_d)
    
    # 记录反力-位移曲线用的数组
    reaction_forces = Float64[]
    displacements = Float64[]
    
    # 建立 VTK 输出文件夹
    mkpath("data/sims")
    
    # --- 3. 开启增量载荷步循环 ---
    println("开始 Staggered 交错求解，总步数: $n_steps")
    for step in 1:n_steps
        # 更新当前的位移加载幅值
        current_disp = (step / n_steps) * max_u_disp
        update!(ch_u, step / n_steps) 
        
        # 将上一载荷步的解作为本步预测值
        u_n .= u_prev
        d_n .= d_prev
        apply!(u_n, ch_u)
        apply!(d_n, ch_d)
        
        println("=== 载荷步 $step / $n_steps | 顶方位移: $(round(current_disp, digits=5)) ===")
        
        # --- 4. 开启内层 Staggered 交错循环 ---
        for iter in 1:max_iter
            # ==========================================
            # 步骤 A: 固定 d, 求解非线性位移场 u
            # ==========================================
            # 牛顿迭代求解 u
            newton_iter = 0
            u_residual_norm = 1.0
            
            while u_residual_norm > tol && newton_iter < 10
                assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
                # 应用 u 的位移边界条件 (直接修改矩阵和残差)
                apply_zero!(K_u, R_u, ch_u)
                
                # R_u 是内力，方程应为 K_u * Δu = -R_u 
                # 但根据 assemble_u 的定义我们可能需要检查正负号
                # 简单起见，标准的 Newton 形式为 K * Δu = -R
                Δu = K_u \ (-R_u)
                apply_zero!(Δu, ch_u) # 边界节点增量为0
                
                u_n .+= Δu
                
                # 重新计算残差的范数作为收敛判断
                assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
                apply_zero!(K_u, R_u, ch_u)
                u_residual_norm = norm(R_u)
                
                newton_iter += 1
            end
            
            # ==========================================
            # 步骤 B: 更新不可逆历史变量 H
            # ==========================================
            update_history!(H_history, dh_u, u_n, mat, cv_u)
            
            # ==========================================
            # 步骤 C: 固定 u (其实是根据 H)，求解线性相场 d
            # ==========================================
            d_old_iter = copy(d_n) # 保存本 iter 开始时的相场，用于判断交错是否收敛
            
            assemble_d!(K_d, F_d, dh_d, H_history, mat, cv_d)
            apply!(K_d, F_d, ch_d) # 应用相场预制裂纹的强制边界条件
            
            d_trial = K_d \ F_d
            d_n .= clamp.(max.(d_trial, d_old_iter), 0.0, 1.0)
            apply!(d_n, ch_d)
            
            # ==========================================
            # 步骤 D: 检查交错循环的收敛性
            # ==========================================
            # 判断准则：这一轮更新的相场，与上一轮相比变化极小
            d_error = norm(d_n - d_old_iter)
            
            println("  - 交错迭代 $iter: Newton_iters=$newton_iter, Δd_error=$(round(d_error, sigdigits=4))")
            
            if d_error < tol
                break # 交错循环收敛，跳出内层 for
            end
            if iter == max_iter
                @warn "  载荷步 $step 在最大交错迭代次数内未收敛！"
            end
        end # 内层交错循环结束
        
        # 确认收敛，更新保存上一载荷步的变量
        u_prev .= u_n
        d_prev .= d_n
        
        # ==========================================
        # 步骤 E: 计算反力与输出 VTK
        # ==========================================
        # 计算整张网格的内力 (包含所有未被零化边界条件破坏的力)
        assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
        
        # 提取顶部竖向位移自由度的反力并求和，避免把底边支座反力混入曲线。
        f_y_total = sum(R_u[dof] for dof in top_y_dofs)
        
        push!(displacements, current_disp)
        push!(reaction_forces, f_y_total)
        
        # 每隔 5 步输出一次 VTK 以节约硬盘
        if step % 5 == 0 || step == n_steps
            VTKGridFile("data/sims/fracture_step_$step", dh_u) do vtk
                write_solution(vtk, dh_u, u_n)
                write_solution(vtk, dh_d, d_n)
            end
        end
        
    end # 外层载荷增量循环结束
    
    println("仿真结束！VTK 文件保存在 data/sims/ 目录下。")
    return displacements, reaction_forces
end
