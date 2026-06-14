"""
    各种各样的实用工具函数，用于简化代码，提升代码可读性和维护性。
"""

"""
    计算 g(d, d_prev, M_d, ρ) = || d - d_prev ||_M^2 - ρ^2 <= 0
"""
function compute_g(d::Vector{Float64}, d_prev::Vector{Float64}, M_d::SparseMatrixCSC{Float64, Int}, ρ::Float64)
    diff_d = d .- d_prev
    M_diff = M_d * diff_d
    g = dot(diff_d, M_diff) - ρ^2
    return g
end

"""
    update_history!(H, dh_u, u_global, mat, cv_u)

根据当前的位移场，更新每个积分点上的相场不可逆历史变量 H。
对应论文公式 (18)。
"""
function update_history!(
    H::Vector{Float64}, 
    dh_u::DofHandler, u_global::Vector{Float64}, 
    mat::PhaseFieldMaterial, cv_u::CellValues
)
    qp_count = 1 # 全局积分点计数器
    for cell in CellIterator(dh_u)
        reinit!(cv_u, cell)
        u_loc = u_global[celldofs(cell)]
        for q_point in 1:getnquadpoints(cv_u)
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            
            # 仅计算拉伸能量部分
            ψ_plus = tensile_energy_density(ε_q, mat)
            
            # 历史变量是递增的 (不可逆性)
            H[qp_count] = max(H[qp_count], ψ_plus)
            
            qp_count += 1
        end
    end
end

"""
计算相场演化的驱动力。
如果 enforce_irreversibility = true，则取历史最大值 (L型试件)；
如果 enforce_irreversibility = false，则直接使用当前应变能 (CT型试件)。
"""
function compute_driving_force!(
    driving_force::Vector{Float64}, 
    dh_u::DofHandler, u_global::Vector{Float64}, 
    mat::PhaseFieldMaterial, cv_u::CellValues,
    enforce_irreversibility::Bool
)
    qp_count = 1
    for cell in CellIterator(dh_u)
        reinit!(cv_u, cell)
        u_loc = u_global[celldofs(cell)]
        
        for qp in 1:getnquadpoints(cv_u)
            # 计算当前积分点的应变
            ε_q = function_symmetric_gradient(cv_u, qp, u_loc)
            
            # 调用你写的拉伸应变能密度函数
            Ψ_plus = tensile_energy_density(ε_q, mat)
            
            # 核心分支判断：
            if enforce_irreversibility
                # L型试件：取历史最大值
                driving_force[qp_count] = max(driving_force[qp_count], Ψ_plus)
            else
                # CT型试件：当前即为驱动力，允许减小
                driving_force[qp_count] = Ψ_plus
            end
            
            qp_count += 1
        end
    end
end

"""
    get_right_dofs(grid, dh_u, dir; tol=1e-12)
提取位于右边界的节点对应的特定方向位移自由度编号，用于计算反力。
dir = 1 代表水平方向(x), dir = 2 代表竖向(y)。
"""
function get_right_dofs(grid, dh_u, dir::Int; tol=1e-12)
    @assert dir == 1 || dir == 2 "dir 必须是 1 (x方向) 或 2 (y方向)"
    
    node_dofs_u = zeros(Int, 2, getnnodes(grid))
    for cell_id in 1:getncells(grid)
        cell = getcells(grid, cell_id)
        dofs = celldofs(dh_u, cell_id)
        for (local_node, node_id) in pairs(cell.nodes)
            node_dofs_u[1, node_id] = dofs[(local_node - 1) * 2 + 1]
            node_dofs_u[2, node_id] = dofs[(local_node - 1) * 2 + 2]
        end
    end

    # 找到最右侧边界的 x 坐标
    coords_x = [node.x[1] for node in grid.nodes]
    right_x = maximum(coords_x)
    
    # 筛选出位于右边界的节点
    right_nodes = findall(x -> isapprox(x, right_x; atol=tol), coords_x)
    
    # 根据传入的 dir 参数提取对应的自由度
    right_dofs = [node_dofs_u[dir, node_id] for node_id in right_nodes]
    
    return unique(right_dofs)
end

"""
    compute_reaction_forces()

计算F_{reaction}, 用于提取反力。
"""
function compute_reaction_forces(f_reac_dof, K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
    # 计算整张网格的内力 (包含所有未被零化边界条件破坏的力)
    assemble_u!(K_u, R_u, dh_u, dh_d, u_n, d_n, mat, cv_u, cv_d)
        
    # 提取右边竖向位移自由度的反力并求和，得到总反力
    f_reac = sum(R_u[dof] for dof in f_reac_dof)
    return f_reac
end

"""
专为 Monolithic 设计的历史变量更新函数
"""
function update_history_mono!(
    H::Vector{Float64}, dh::DofHandler, x_global::Vector{Float64}, 
    mat::PhaseFieldMaterial, cv_u::CellValues
)
    qp_count = 1
    u_range = dof_range(dh, :u)
    for cell in CellIterator(dh)
        reinit!(cv_u, cell)
        u_loc = x_global[celldofs(cell)][u_range]
        for q_point in 1:getnquadpoints(cv_u)
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            ψ_plus = tensile_energy_density(ε_q, mat)
            H[qp_count] = max(H[qp_count], ψ_plus)
            qp_count += 1
        end
    end
end
function compute_driving_force_mono!(
    driving_force::Vector{Float64}, 
    dh::DofHandler, x_global::Vector{Float64}, 
    mat::PhaseFieldMaterial, cv_u::CellValues,
    enforce_irreversibility::Bool
)
    qp_count = 1
    u_range = dof_range(dh, :u)
    for cell in CellIterator(dh)
        reinit!(cv_u, cell)
        u_loc = x_global[celldofs(cell)][u_range]
        
        for qp in 1:getnquadpoints(cv_u)
            # 计算当前积分点的应变
            ε_q = function_symmetric_gradient(cv_u, qp, u_loc)

            Ψ_plus = tensile_energy_density(ε_q, mat)
            
            # 核心分支判断：
            if enforce_irreversibility
                # L型试件：取历史最大值
                driving_force[qp_count] = max(driving_force[qp_count], Ψ_plus)
            else
                # CT型试件：当前即为驱动力，允许减小
                driving_force[qp_count] = Ψ_plus
            end 
            qp_count += 1
        end
    end
end

"""
计算 Γ 约束残差 f_Γ = 𝒢_f(d_curr) - 𝒢_f(d_prev) - ρ
以及约束方程关于全体自由度的梯度向量 K_λa。
"""
function evaluate_gamma_constraint(
    dh::DofHandler, a_global::Vector{Float64}, G_prev::Float64,
    mat::PhaseFieldMaterial, cv_d::CellValues, ρ::Float64
)
    K_λa = zeros(Float64, ndofs(dh))
    G_curr = 0.0
    
    n_basefuncs_d = getnbasefunctions(cv_d)
    K_λd_loc = zeros(Float64, n_basefuncs_d)
    
    d_range = dof_range(dh, :d)

    for cell in CellIterator(dh)
        reinit!(cv_d, cell)
        global_dofs = celldofs(cell)
        
        a_loc = a_global[global_dofs]
        d_loc = a_loc[d_range]
        
        fill!(K_λd_loc, 0.0)
        
        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            d_q = function_value(cv_d, q_point, d_loc)
            ∇d_q = function_gradient(cv_d, q_point, d_loc)
            
            # 累加表面能
            G_curr += (mat.gc / (2 * mat.l)) * (d_q^2 + mat.l^2 * (∇d_q ⋅ ∇d_q)) * dΩ
            
            # 计算对 d 的局部梯度
            for i in 1:n_basefuncs_d
                N_i = shape_value(cv_d, q_point, i)
                ∇N_i = shape_gradient(cv_d, q_point, i)
                K_λd_loc[i] += (mat.gc / mat.l) * (d_q * N_i + mat.l^2 * (∇d_q ⋅ ∇N_i)) * dΩ
            end
        end
        
        # 组装到全局向量
        for i in 1:n_basefuncs_d
            global_idx = global_dofs[d_range[i]]
            K_λa[global_idx] += K_λd_loc[i]
        end
    end
    
    f_Γ = G_curr - G_prev - ρ
    return f_Γ, K_λa, G_curr
end

"""
装配 H1 弧长法专用的恒定几何矩阵 H = M + l^2 S
"""
function assemble_H1_matrix(dh::DofHandler, cv_d::CellValues, l::Float64)
    # 使用 allocate_matrix 替代已弃用的 create_sparsity_pattern
    K_sparse = allocate_matrix(dh)
    assembler = start_assemble(K_sparse)
    
    # 提取相场 d 的局部自由度范围
    d_range = dof_range(dh, :d)
    n_d = length(d_range)

    for cell in CellIterator(dh)
        reinit!(cv_d, cell)
        # 尺寸与该单元的所有自由度（u和d）相同，但只填充 d 的部分
        Ke = zeros(ndofs_per_cell(dh), ndofs_per_cell(dh))

        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            for i in 1:n_d
                N_i = shape_value(cv_d, q_point, i)
                ∇N_i = shape_gradient(cv_d, q_point, i)
                for j in 1:n_d
                    N_j = shape_value(cv_d, q_point, j)
                    ∇N_j = shape_gradient(cv_d, q_point, j)
                    
                    # H1 范数的核： N_i * N_j + l^2 * (∇N_i ⋅ ∇N_j)
                    # 放入 d_range 对应的行列位置
                    Ke[d_range[i], d_range[j]] += (N_i * N_j + l^2 * (∇N_i ⋅ ∇N_j)) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Ke)
    end
    return K_sparse
end

"""
装配 L2 弧长法专用的恒定几何矩阵 M (相当于相场的质量矩阵)
"""
function assemble_L2_matrix(dh::DofHandler, cv_d::CellValues)
    K_sparse = allocate_matrix(dh)
    assembler = start_assemble(K_sparse)
    
    d_range = dof_range(dh, :d)
    n_d = length(d_range)

    for cell in CellIterator(dh)
        reinit!(cv_d, cell)
        Ke = zeros(ndofs_per_cell(dh), ndofs_per_cell(dh))

        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            for i in 1:n_d
                N_i = shape_value(cv_d, q_point, i)
                # (不需要求 ∇N_i 了)
                for j in 1:n_d
                    N_j = shape_value(cv_d, q_point, j)
                    
                    # L2 范数的核：纯纯的 N_i * N_j
                    Ke[d_range[i], d_range[j]] += (N_i * N_j) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Ke)
    end
    return K_sparse
end


"""
M-EM 求解器的不同策略集下的迭代求解函数
"""
# ====================================================================
# [内部函数] 无约束 Newton (Eq 47)
# ====================================================================
function solve_newton_inactive!(
    a_cur, K_mono, r_mono, dh, ch_zero, H_old, mat, cv_u, cv_d, max_iter, tol
)
    for iter in 1:max_iter
        assemble_monolithic!(K_mono, r_mono, dh, a_cur, H_old, mat, cv_u, cv_d)
        
        r_check = -copy(r_mono)
        apply_zero!(r_check, ch_zero)
        res_norm = norm(r_check)
        
        if res_norm < tol
            return true, a_cur, iter
        end

        apply_zero!(K_mono, r_mono, ch_zero)
        Δa = K_mono \ (-r_mono)
        a_cur .+= Δa
    end
    return false, a_cur, max_iter
end


# ====================================================================
# [内部函数] 增广系统 Newton (Eq 44, 45) - 舒尔补极致优化版
# ====================================================================
function solve_newton_active!(
    a_cur, a_prev, λ_rho, K_mono, r_mono, dh, ch_zero, ch_a, H_old, mat, cv_u, cv_d, 
    M_d, ρ, idx_u, idx_d, max_iter, tol
)
    n_dofs = length(a_cur)
    pdofs = ch_zero.prescribed_dofs
    
    for iter in 1:max_iter
        # 1. 组装标准单块刚度矩阵和残差
        assemble_monolithic!(K_mono, r_mono, dh, a_cur, H_old, mat, cv_u, cv_d)
        
        # 2. 计算约束 g_val = (d - d_prev)^T * M_d * (d - d_prev) - ρ^2
        Δd = a_cur[idx_d] .- a_prev[idx_d]
        M_Δd = M_d * Δd
        g_val = dot(Δd, M_Δd) - ρ^2
        
        # 约束梯度 dg_da (在全局自由度下)
        dg_da = zeros(n_dofs)
        dg_da[idx_d] .= 2.0 .* M_Δd
        
        # 3. 构造修正后的标准单块稀疏矩阵 K_aug_base = K_mono + [0 0; 0 2*λ_rho*M_d]
        K_aug_base = copy(K_mono)
        I, J, V = findnz(M_d)
        for k in eachindex(V)
            # 仅在非零元素上叠加修改，保持稀疏性
            K_aug_base[idx_d[I[k]], idx_d[J[k]]] += 2.0 * λ_rho * V[k]
        end
        
        # 4. 检查收敛性 (计算增广系统残差)
        r_aug_a = r_mono .+ λ_rho .* dg_da
        
        r_check = zeros(n_dofs + 1)
        r_check[1:n_dofs] .= r_aug_a
        r_check[end] = g_val
        r_check[pdofs] .= 0.0  # 忽略 Dirichlet 约束自由度
        
        res_norm = norm(r_check)
        if res_norm < tol
            return true, a_cur, λ_rho, iter
        end

        # 5. 施加齐次 Dirichlet 边界条件 (利用 Ferrite 内置的高效方法)
        # 这会将 K_aug_base 中 pdofs 对应的行和列置零，对角线置 1
        apply_zero!(K_aug_base, r_aug_a, ch_zero)
        
        # 约束梯度也必须在 Dirichlet 自由度处清零，以保持线性系统的相容性
        dg_da_bc = copy(dg_da)
        dg_da_bc[pdofs] .= 0.0

        # 6. 【核心】通过舒尔补求解，不破坏大矩阵的稀疏结构
        # 求解两个 N*N 的稀疏线性方程组
        v_r = K_aug_base \ (-r_aug_a)
        v_g = K_aug_base \ (-dg_da_bc)
        
        # 求解对偶变量增量 Δλ
        denom = dot(dg_da_bc, v_g)
        if abs(denom) < 1e-12
            denom = sign(denom) * 1e-12 # 防止分母为 0 导致数值崩溃
        end
        Δλ = (-g_val - dot(dg_da_bc, v_r)) / denom
        
        # 计算位移/相场增量 (Dirichlet 自由度处的增量自动为 0)
        Δa = v_r .+ Δλ .* v_g
        
        # 7. 更新变量
        a_cur .+= Δa
        λ_rho += Δλ
    end
    
    return false, a_cur, λ_rho, max_iter
end