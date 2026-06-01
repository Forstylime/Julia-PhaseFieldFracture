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
连续收敛失败时减半 ρ；每次成功后逐步恢复（×1.1，上限为初始值）。
"""
function adapt_rho!(ρ::Float64, ρ_init::Float64, success::Bool)
    if success
        return min(ρ * 1.2, ρ_init)
    else
        return ρ * 0.5
    end
end