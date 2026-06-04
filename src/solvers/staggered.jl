# src/solvers/staggered.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
    solve_staggered(setup::TensionSetup, mat::PhaseFieldMaterial;
                    n_steps = 100, tol = 1e-4)

执行标准交错求解法 (Staggered Scheme)。
在每个载荷步中，交错求解位移场和相场，直到残差收敛。
位移加载幅值从 setup.final_displacement 读取。
"""
function solve_staggered(
    setup::TensionSetup, mat::PhaseFieldMaterial;
    n_steps = 100,          # 总载荷步数 (论文设为 100 步)
    tol = 1e-5,             # Newton 迭代和交错循环的收敛容差
    max_iter = 20           # 内层交错循环的最大允许次数
)
    # --- 1. 提取网格与自由度 ---
    grid = setup.grid
    dh_u = setup.dh_u
    dh_d = setup.dh_d
    ch_u = setup.ch_u
    ch_d = setup.ch_d

    right_x_dofs = get_right_dofs(grid, dh_u)
    
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
    
    # 记录反力-位移曲线和能量演化用的数组
    reaction_forces = Float64[0]
    displacements = Float64[0]
    elastic_energies = Float64[0]
    surface_energies = Float64[0]
    
    # 建立 VTK 输出文件夹
    mkpath("data/sims")

    # 计时和迭代计数
    t_start = time()
    total_newton_iters = 0

    # --- 3. 开启增量载荷步循环 ---
    println("开始 Staggered 交错求解，总步数: $n_steps")
    for step in 1:n_steps
        # 更新当前的位移加载幅值
        current_disp = (step / n_steps) * setup.final_displacement
        update!(ch_u, step / n_steps) 
        
        # 将上一载荷步的解作为本步预测值
        u_n .= u_prev
        d_n .= d_prev
        apply!(u_n, ch_u)
        apply!(d_n, ch_d)
        
        println("=== 载荷步 $step / $n_steps | 位移: $(round(current_disp, digits=5)) ===")
        
        # --- 4. 开启内层 Staggered 交错循环 ---
        for iter in 1:max_iter
            # ==========================================
            # 步骤 A: 固定 d, 求解非线性位移场 u
            # ==========================================
            # 牛顿迭代求解 u
            newton_iter = 0
            u_residual_norm = 1.0
            
            # 1. 进入循环前，先进行首次组装（假设此时 u_n 已应用了当前载荷步的非零边界条件）
            assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
            apply_zero!(K_u, R_u, ch_u)
            u_residual_norm = norm(R_u)

            while u_residual_norm > tol && newton_iter < 10
                # 2. 求解位移增量（注意：-R_u 会产生临时分配，若追求极致性能可用 ldiv! 并提前取反）
                Δu = K_u \ (-R_u)
                apply_zero!(Δu, ch_u) # 确保边界处的增量严格为 0
                
                # 3. 更新位移
                u_n .+= Δu
                
                # 4. 在新的位移下组装，为下一次迭代（或退出条件判断）做准备
                assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
                apply_zero!(K_u, R_u, ch_u)
                
                # 5. 更新残差范数和迭代步数
                u_residual_norm = norm(R_u)
                newton_iter += 1
            end
            total_newton_iters += newton_iter
            
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
            
            if mod(iter, 5) == 0 
                @info " - 交错迭代 $iter, Δd_error=$(round(d_error, sigdigits=4))"
            end
            
            if d_error < tol
                break # 交错循环收敛，跳出内层 for
            end
            if iter == max_iter
                @warn " - 载荷步 $step 在最大交错迭代次数内未收敛！- "
            end
        end # 内层交错循环结束
        
        # 确认收敛，更新保存上一载荷步的变量
        u_prev .= u_n
        d_prev .= d_n

        # 计算当前步的弹性体能量和断裂表面能
        push!(elastic_energies, elastic_energy(dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d))
        push!(surface_energies, surface_energy(dh_d, d_n, mat, cv_d))

        # ==========================================
        # 步骤 E: 计算反力与输出 VTK
        # ==========================================
        f__total = compute_reaction_forces(right_x_dofs, K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
        
        push!(displacements, current_disp)
        push!(reaction_forces, f__total)
        
        # 每隔 5 步输出一次 VTK 以节约硬盘
        if step % 5 == 0 || step == n_steps
            VTKGridFile("data/sims/staggered/fracture_step_$step", dh_u) do vtk
                write_solution(vtk, dh_u, u_n)
                write_solution(vtk, dh_d, d_n)
            end
        end
        
    end # 外层载荷增量循环结束
    
    println("仿真结束！VTK 文件保存在 data/sims/staggered 目录下。")
    println("总Newton迭代次数: $total_newton_iters")
    println("计算耗时: $(round(time() - t_start, digits=2)) 秒")
    return displacements, reaction_forces, elastic_energies, surface_energies
end
