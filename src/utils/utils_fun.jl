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