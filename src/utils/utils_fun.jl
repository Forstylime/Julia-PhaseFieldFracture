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
    get_right_dofs(grid, dh_u, tol)
提取位于右边界的节点对应的竖向位移自由度编号，用于计算反力。
"""
function get_right_dofs(grid, dh_u, tol=1e-12)
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
    right_nodes = findall(x -> isapprox(x, right_x; atol=tol), coords_x)
    right_dofs = [node_dofs_u[2, node_id] for node_id in right_nodes]
    
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

"""
    adapt_rho!(ρ, success)

Crisfield 弧长半径自适应。
连续收敛失败时减半 ρ；每次成功后逐步恢复（×1.2，上限为初始值）。
"""
function adapt_rho!(ρ::Float64, ρ_init::Float64, success::Bool)
    if success
        return min(ρ * 1.2, ρ_init)
    else
        return ρ * 0.5
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