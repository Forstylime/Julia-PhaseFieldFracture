# src/physics/energies.jl

using Ferrite

"""
    elastic_energy(dh_u, dh_d, u_global, d_global, mat, cv_u, cv_d)

计算弹性体能量 Ψ(u, d) = ∫_Ω [g(d) ψ₀⁺(ε) + ψ₀⁻(ε)] dΩ。

在位移-控制加载条件下，总势能 ℱ = Ψ + 𝒢_f（无外力功）。
"""
function elastic_energy(
    dh_u::DofHandler, dh_d::DofHandler,
    u_global::Vector{Float64}, d_global::Vector{Float64},
    mat::PhaseFieldMaterial,
    cv_u::CellValues, cv_d::CellValues,
)
    energy = 0.0
    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)
        u_loc = u_global[celldofs(cell_u)]
        d_loc = d_global[celldofs(cell_d)]
        for q_point in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, q_point)
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            d_q = function_value(cv_d, q_point, d_loc)
            energy += elastic_energy_density(ε_q, d_q, mat) * dΩ
        end
    end
    return energy
end

"""
    surface_energy(dh_d, d_global, mat, cv_d)

计算断裂表面能 𝒢_f(d) = ∫_Ω g_c/(2l) [d² + l²|∇d|²] dΩ。
对应论文 Eq. (5)。
"""
function surface_energy(
    dh_d::DofHandler, d_global::Vector{Float64},
    mat::PhaseFieldMaterial, cv_d::CellValues,
)
    energy = 0.0
    for cell in CellIterator(dh_d)
        reinit!(cv_d, cell)
        d_loc = d_global[celldofs(cell)]
        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            d_q = function_value(cv_d, q_point, d_loc)
            ∇d_q = function_gradient(cv_d, q_point, d_loc)
            energy += (mat.gc / (2 * mat.l)) * (d_q^2 + mat.l^2 * (∇d_q ⋅ ∇d_q)) * dΩ
        end
    end
    return energy
end

"""
    total_energy(dh_u, dh_d, u_global, d_global, mat, cv_u, cv_d)

计算系统总能量 ℱ(u, d) = Ψ(u, d) + 𝒢_f(d)。

当前算例采用位移-控制加载（Dirichlet BC），无外力功 Π_ext = 0。
"""
function total_energy(
    dh_u::DofHandler, dh_d::DofHandler,
    u_global::Vector{Float64}, d_global::Vector{Float64},
    mat::PhaseFieldMaterial,
    cv_u::CellValues, cv_d::CellValues,
)
    return elastic_energy(dh_u, dh_d, u_global, d_global, mat, cv_u, cv_d) +
           surface_energy(dh_d, d_global, mat, cv_d)
end



"""
    elastic_energy_monolithic(dh::DofHandler, x_global::Vector{Float64}, mat, cv_u, cv_d)

整体法 (Monolithic) 的弹性能计算版本。
由于 u 和 d 存储在同一个状态向量 x_global 中，我们需要使用 dof_range 提取局部向量。
"""
function elastic_energy_monolithic(
    dh::DofHandler, x_global::Vector{Float64},
    mat::PhaseFieldMaterial,
    cv_u::CellValues, cv_d::CellValues,
)
    energy = 0.0
    # 获取各个场在单元内的局部自由度索引范围
    u_range = Ferrite.dof_range(dh, :u)
    d_range = Ferrite.dof_range(dh, :d)

    for cell in CellIterator(dh)
        reinit!(cv_u, cell)
        reinit!(cv_d, cell)
        
        # 提取当前单元包含所有场（u和d）的局部向量
        x_loc = x_global[celldofs(cell)]
        
        # 拆分出位移和相场的局部向量
        u_loc = x_loc[u_range]
        d_loc = x_loc[d_range]
        
        for q_point in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, q_point)
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            d_q = function_value(cv_d, q_point, d_loc)
            energy += elastic_energy_density(ε_q, d_q, mat) * dΩ
        end
    end
    return energy
end

"""
    surface_energy_monolithic(dh::DofHandler, x_global::Vector{Float64}, mat, cv_d)

整体法 (Monolithic) 的断裂表面能计算版本。
"""
function surface_energy_monolithic(
    dh::DofHandler, x_global::Vector{Float64},
    mat::PhaseFieldMaterial, cv_d::CellValues,
)
    energy = 0.0
    d_range = Ferrite.dof_range(dh, :d)

    for cell in CellIterator(dh)
        reinit!(cv_d, cell)
        x_loc = x_global[celldofs(cell)]
        d_loc = x_loc[d_range]
        
        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            d_q = function_value(cv_d, q_point, d_loc)
            ∇d_q = function_gradient(cv_d, q_point, d_loc)
            energy += (mat.gc / (2 * mat.l)) * (d_q^2 + mat.l^2 * (∇d_q ⋅ ∇d_q)) * dΩ
        end
    end
    return energy
end

"""
    total_energy_monolithic(dh, x_global, mat, cv_u, cv_d)
"""
function total_energy_monolithic(
    dh::DofHandler, x_global::Vector{Float64},
    mat::PhaseFieldMaterial, cv_u::CellValues, cv_d::CellValues,
)
    return elastic_energy_monolithic(dh, x_global, mat, cv_u, cv_d) +
           surface_energy_monolithic(dh, x_global, mat, cv_d)
end