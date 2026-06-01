# src/solvers/crisfield.jl

using Ferrite
using LinearAlgebra
using SparseArrays

"""
统一 DofHandler 架构下的边界条件施加函数（处理非齐次/齐次解析）。
"""
function apply_arc_length_bc!(K::SparseMatrixCSC{Float64, Int}, f::Vector{Float64}, 
                              ch::ConstraintHandler, a_bc::Vector{Float64}, apply_zero::Bool)
    n_dofs = size(K, 1)
    is_constrained = falses(n_dofs)
    is_constrained[ch.prescribed_dofs] .= true
    
    if !apply_zero
        # 把非对角项移动到右端项 f
        f .-= K * a_bc
        for i in ch.prescribed_dofs
            f[i] = a_bc[i]
        end
    else
        # 残差方程，边界增量为0
        for i in ch.prescribed_dofs
            f[i] = 0.0
        end
    end
    
    # 极速将约束行列清零，对角线置1
    for j in 1:size(K, 2)
        for k in K.colptr[j]:(K.colptr[j+1]-1)
            i = K.rowval[k]
            if is_constrained[i] || is_constrained[j]
                K.nzval[k] = (i == j) ? 1.0 : 0.0
            end
        end
    end
    
    return K, f
end

"""
Crisfield 弧长法求解器 (Monolithic 架构版)
输入的是 `MonolithicTensionSetup`
"""
function solve_crisfield(
    setup::MonolithicTensionSetup, mat::PhaseFieldMaterial;
    ρ_init::Float64 = 0.05,
    ρ_min::Float64 = 1e-6,    # ρ 低于此值视为收敛停滞
    max_steps::Int = 200,
    tol::Float64 = 1e-6,
    max_newton::Int = 15,
    output_freq::Int = 5,
    λ_max::Float64 = 1.0,
)
    grid = setup.grid
    dh = setup.dh
    ch = setup.ch
    n_total = ndofs(dh)

    # ================================================================
    # 获取 u 和 d 在全局向量中的索引掩码 (因为它们交织在一起)
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

    # 提取右端点位移自由度（算反力用）
    coords_x = [node.x[1] for node in grid.nodes]
    right_x = maximum(coords_x)
    right_dofs = Int[]
    for cell in CellIterator(dh)
        for (i, node_id) in enumerate(cell.nodes)
            if isapprox(grid.nodes[node_id].x[1], right_x; atol = 1e-8)
                push!(right_dofs, celldofs(cell)[u_range[(i-1)*2 + 2]]) # 获取 uy (竖向加载)
            end
        end
    end
    right_dofs = unique(right_dofs)

    # ================================================================
    # 构建 λ = 1.0 时的全局参考向量 a_bc_ref
    # ================================================================
    a_bc_ref = zeros(n_total)
    update!(ch, 1.0)
    apply!(a_bc_ref, ch)
    update!(ch, 0.0) # 恢复

    # 初始化状态 (统一向量)
    x_n = zeros(n_total)
    x_cur = zeros(n_total)
    x_prev = zeros(n_total) # 记录上一步的收敛状态
    λ_n = 0.0; λ_cur = 0.0
    λ_max_reached = 0.0    # 跟踪最大 λ (检测振荡用)

    qr = QuadratureRule{RefQuadrilateral}(2)
    cv_u = CellValues(qr, Lagrange{RefQuadrilateral, 1}()^2)
    cv_d = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    n_qpoints = getncells(grid) * getnquadpoints(cv_u)
    H_old = zeros(n_qpoints)

    # 只在这里分配一次巨大的稀疏矩阵拓扑结构！
    K_mono = allocate_matrix(dh)
    r_mono = zeros(n_total)
    
    displacements = Float64[0.0]
    reaction_forces = Float64[0.0]
    elastic_energies = Float64[0.0]
    surface_energies = Float64[0.0]

    num = 0.0
    den = 0.0
    δλ = 0.0

    mkpath("data/sims/crisfield")
    ρ = ρ_init
    n_success = 0

    println("开始 Crisfield 弧长法 (Monolithic)，ρ_init = $ρ, λ_max = $λ_max")

    while λ_n < λ_max && n_success < max_steps
        println("=== 弧长步 $n_success | λ = $(round(λ_n, digits=6)) | ρ = $(round(ρ, digits=6)) ===")

        # =============================================
        # PHASE 1: PREDICTOR (预测步)
        # =============================================
        assemble_monolithic!(K_mono, r_mono, dh, x_n, H_old, mat, cv_u, cv_d)

        K_f_pred = copy(K_mono)
        f_f_pred = zeros(n_total)
        apply_arc_length_bc!(K_f_pred, f_f_pred, ch, a_bc_ref, false)
        
        δa_f_pred = K_f_pred \ f_f_pred
        
        Nu = length(idx_u)
        δu_f_pred = δa_f_pred[idx_u] 
        norm_f_scaled = norm(δu_f_pred) / sqrt(Nu)

        if norm_f_scaled < 1e-14
            @warn "刚度退化为0，计算终止。"
            break
        end

        # 预测方向：基于前一步全状态增量 (u + d), 避免 peak 附近 Δu≈0 导致符号翻转
        if n_success == 0
            sign_λ = 1.0
        else
            Δa_prev = x_n .- x_prev
            dot_prev = dot(Δa_prev, δa_f_pred)
            sign_λ = abs(dot_prev) < 1e-16 ? 1.0 : sign(dot_prev)
        end
        
        Δλ_pred = sign_λ * ρ / norm_f_scaled
        λ_cur = λ_n + Δλ_pred

        x_cur .= x_n .+ Δλ_pred .* δa_f_pred

        # 同步边界解析值
        for i in ch.prescribed_dofs
            x_cur[i] = λ_cur * a_bc_ref[i]
        end

        # =============================================
        # PHASE 2: CORRECTOR (Monolithic Newton, 论文公式 24-26)
        # =============================================
        converged = false
        n_newton = 0
        for iter in 1:max_newton
            n_newton = iter
            assemble_monolithic!(K_mono, r_mono, dh, x_cur, H_old, mat, cv_u, cv_d)

            # 1. 计算位移增量 Δu_k 和柱面弧长约束残差 (标准二次型 Crisfield)
            Δu_k = x_cur[idx_u] .- x_n[idx_u]
            f_C = dot(Δu_k, Δu_k) / Nu - ρ^2   # ψ = ‖Δu‖²/Nu - ρ²

            # 收敛检查：消除位移边界反力后的力残差 + 弧长约束残差
            r_check = copy(r_mono)
            r_check[ch.prescribed_dofs] .= 0.0

            r_norm_u = norm(r_check[idx_u]) / sqrt(Nu)
            r_norm_d = norm(r_check[idx_d]) / sqrt(length(idx_d))

            if r_norm_u < tol && r_norm_d < tol && abs(f_C) < tol
                converged = true
                break
            end

            # 2. 载荷方向增量 δa_f
            K_f = copy(K_mono)
            f_f = zeros(n_total)
            apply_arc_length_bc!(K_f, f_f, ch, a_bc_ref, false)
            δa_f = K_f \ f_f
            δu_f = δa_f[idx_u]

            # 3. 残差方向增量 δa_r
            K_r = copy(K_mono)
            f_r = -copy(r_mono)
            apply_arc_length_bc!(K_r, f_r, ch, a_bc_ref, true)
            δa_r = K_r \ f_r
            δu_r = δa_r[idx_u]

            # 4. 根据论文公式 (26) 求解 Δλ
            # 约束梯度: ∇ψ = 2·Δu_k / Nu  (二次型 ψ = ‖Δu‖²/Nu - ρ² 的梯度)
            grad_ψ = 2.0 .* Δu_k ./ Nu

            # ψ + ∇ψ·δa_r = f_C + grad_ψ·δu_r
            num = f_C + dot(grad_ψ, δu_r)
            # ∇ψ·δa_f + ∂ψ/∂λ (= 0 for cylindrical Crisfield)
            den = dot(grad_ψ, δu_f)

            if abs(den) < 1e-14
                @warn "Crisfield: singular denominator (∇ψ · δu_f ≈ 0), aborting step"
                converged = false
                break
            end
            δλ = -num / den

            # 5. 更新状态 (论文公式 25 第一项)
            x_cur .+= δa_r .+ δλ .* δa_f
            λ_cur += δλ

            # 6. 强制 d 不可逆，可能不需要，暂时注释
            #for i in idx_d
            #    x_cur[i] = clamp(x_cur[i], x_n[i], 1.0)
            #end

            # 7. 强制拉回边界约束
            for i in ch.prescribed_dofs
                x_cur[i] = λ_cur * a_bc_ref[i]
            end
        end

        # =============================================
        # PHASE 3: 收敛后处理 → 步长控制
        # =============================================
        if !converged
            if ρ <= ρ_min
                @warn "ρ 已降至最小值 ($(ρ_min))，无法继续收敛。终止求解。"
                break
            end
            @warn "未在 $max_newton 次迭代内收敛。回退并减半弧长 ρ。"
            ρ *= 0.5
            continue
        end

        # 收敛后强制 d 不可逆性 (安全兜底, 已在 Newton 迭代中执行)
        for i in idx_d
            x_cur[i] = clamp(x_cur[i], x_n[i], 1.0)
        end

        n_success += 1

        # 基于 Newton 迭代次数自适应步长:
        # n_newton ≤ 4  → ρ 放大 (步长可增大)
        # n_newton ≥ 8  → ρ 缩小 (步长过大, 接近收敛极限)
        # 4 < n_newton < 8 → ρ 保持不变
        if n_newton <= 4
            ρ = min(ρ * 1.5, ρ_init * 3.0)
        elseif n_newton >= 8
            ρ = max(ρ * 0.7, ρ_min * 10.0)
        end

        # 更新历史变量
        update_history_mono!(H_old, dh, x_cur, mat, cv_u)

        x_prev .= x_n
        x_n .= x_cur
        λ_n = λ_cur

        # 振荡检测: 若 d 已饱和且 λ 远低于历史最大值 → 断裂已完全发展, 结束
        λ_max_reached = max(λ_max_reached, λ_n)
        d_max = maximum(x_n[idx_d])
        if d_max > 0.95 && λ_n < 0.15 * λ_max_reached
            @info "断裂已完全发展 (d_max=$d_max), λ=$λ_n 远低于峰值 λ=$(round(λ_max_reached, digits=4))。终止求解。"
            break
        end
        @info "max d = $d_max"

        # 计算并保存反力 (直接拿刚更新完的 r_mono 里的值)
        # 注意：上面最后一次循环跳出时并没有计算最新状态的残差，为求精确反力再组装一次
        assemble_monolithic!(K_mono, r_mono, dh, x_n, H_old, mat, cv_u, cv_d)
        f_total = sum(r_mono[dof] for dof in right_dofs)

        push!(displacements, λ_n * setup.final_displacement)
        push!(reaction_forces, f_total)
        push!(elastic_energies, elastic_energy_monolithic(dh, x_n, mat, cv_u, cv_d))
        push!(surface_energies, surface_energy_monolithic(dh, x_n, mat, cv_d)) 

        if n_success % output_freq == 0
            # 使用 Ferrite 原生的 VTK 导出，注意把全局状态 x_n 切分给 u 和 d
            VTKGridFile("data/sims/crisfield/fracture_step_$n_success", setup.grid) do vtk
                write_solution(vtk, dh, x_n) # Ferrite 会自动将 u 和 d 作为不同场写入
            end
        end
    end

    println("Crisfield 仿真结束。VTK 保存在 data/sims/crisfield 目录下。")
    return displacements, reaction_forces, elastic_energies, surface_energies
end