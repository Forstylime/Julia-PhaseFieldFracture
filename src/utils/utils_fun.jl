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
    solve_crisfield_quadratic(a, b, c, Δa_old, δa_r, δa_λ)
Crisfield 二次方程求解器，选择使得新位移增量向量与旧增量向量夹角更小的那个根。
"""
function solve_crisfield_quadratic(a::Float64, b::Float64, c::Float64, 
                                   Δa_old::Vector{Float64}, 
                                   δa_r::Vector{Float64}, 
                                   δa_λ::Vector{Float64})
    
    discriminant = b^2 - 4.0 * a * c
    
    if discriminant < 0
        # 极少数情况下，迭代漂移太远，切线与目标球面无交点。
        # 严谨的做法是抛出异常，让外层程序减小弧长重试。
        error("Crisfield二次方程无实数解 (判别式 < 0)，需要减小弧长或重置本步！")
    end
    
    # 求解两个根
    δλ_1 = (-b + sqrt(discriminant)) / (2.0 * a)
    δλ_2 = (-b - sqrt(discriminant)) / (2.0 * a)
    
    # 根据两个根计算两种可能的“新位移增量向量”
    Δa_new_1 = Δa_old .+ δa_r .+ δλ_1 .* δa_λ
    Δa_new_2 = Δa_old .+ δa_r .+ δλ_2 .* δa_λ
    
    # 计算新向量与旧向量的点积 (cosθ 越大代表夹角越小)
    dot_1 = dot(Δa_old, Δa_new_1)
    dot_2 = dot(Δa_old, Δa_new_2)
    
    # 返回使得点积更大的那个 δλ
    if dot_1 > dot_2
        return δλ_1
    else
        return δλ_2
    end
end