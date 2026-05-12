# src/physics/constitutive.jl

using Tensors
using LinearAlgebra

"""
    PhaseFieldMaterial

存储线弹性相场断裂模型的材料参数。
论文采用 Lamé 常数 λ 和 μ，以及断裂能 gc 和尺度参数 l。
"""
Base.@kwdef struct PhaseFieldMaterial
    E::Float64        = 210e3   # 杨氏模量 (MPa)
    ν::Float64        = 0.3     # 泊松比
    # Lamé 常数 (假设平面应变，若是平面应力后续可在计算能量时修正)
    λ::Float64        = (E * ν) / ((1 + ν) * (1 - 2ν))
    μ::Float64        = E / (2 * (1 + ν))
    # 相场断裂参数
    gc::Float64       = 2.7     # 临界断裂能 (N/mm)
    l::Float64        = 0.01    # 长度尺度参数 (mm)
    k_tol::Float64    = 1e-8    # 残余刚度参数 k0 (防止 d=1 时刚度矩阵奇异)
end

"""
    macauley_plus(x) / macauley_minus(x)

Macaulay 括号，用于提取正数或负数部分，对应论文公式里的 ⟨x⟩±。
"""
macauley_plus(x::Real) = max(x, zero(x))
macauley_minus(x::Real) = min(x, zero(x))

"""
    spectral_decomposition(strain::AbstractMatrix)

兼容矩阵形式应变输入的谱分解接口，返回特征值和特征向量。
该方法保留早期测试和工具函数使用的矩阵 API。
"""
function spectral_decomposition(strain::AbstractMatrix)
    F = eigen(Symmetric(strain))
    return F.values, F.vectors
end

"""
    positive_strain(strain::AbstractMatrix)

返回矩阵应变的正谱部分，用于旧的 `MaterialParameters` 工作流。
"""
function positive_strain(strain::AbstractMatrix)
    values, vectors = spectral_decomposition(strain)
    positive_values = Diagonal(map(v -> max(v, zero(v)), values))
    return vectors * positive_values * vectors'
end

"""
    stress(strain, d, material::MaterialParameters)

基于旧 `MaterialParameters` 类型的拉伸谱分解应力计算。
保留该方法可以避免现有测试和脚本在引入 `PhaseFieldMaterial` 后失效。
"""
function stress(strain::AbstractMatrix, d, material::MaterialParameters)
    eps_pos = positive_strain(strain)
    tr_pos = tr(eps_pos)
    sigma_pos = material.lambda * tr_pos * I + 2 * material.mu * eps_pos
    return residual_degradation(d, material.kappa) * sigma_pos
end

"""
    spectral_decomposition(ε::SymmetricTensor{2, 2})

对 2D 应变张量进行谱分解，返回正部分 ε+ 和负部分 ε-。
对应论文 Eq. (4)。
"""
function spectral_decomposition(ε::SymmetricTensor{2, 2, T}) where T
    # 对对称张量求特征值和特征向量
    # 注意：因为 ForwardDiff 传进来的类型 T 可能是 Dual 数，
    # 这里我们保证类型泛型 T，以支持自动微分穿透！
    eig = eigen(ε)
    λ1, λ2 = eig.values
    v1, v2 = eig.vectors[:, 1], eig.vectors[:, 2]
    
    # 构造特征向量的并矢 (v ⊗ v)
    M1 = symmetric(v1 ⊗ v1)
    M2 = symmetric(v2 ⊗ v2)
    
    # 组装拉伸应变张量 ε+ 和压缩应变张量 ε-
    ε_plus  = macauley_plus(λ1) * M1  + macauley_plus(λ2) * M2
    ε_minus = macauley_minus(λ1) * M1 + macauley_minus(λ2) * M2
    
    return ε_plus, ε_minus
end

"""
    elastic_energy_density(ε::SymmetricTensor{2,2}, d::Real, mat::PhaseFieldMaterial)

计算单个积分点上的纯弹性能密度 Ψ(ε, d)。
结合了微裂纹闭合效应 (MCR-effect)，对应论文 Eq. (2) 和 Eq. (3)。
"""
function elastic_energy_density(ε::SymmetricTensor{2,2,T}, d::Real, mat::PhaseFieldMaterial) where T
    # 1. 谱分解得到拉伸与压缩应变
    ε_plus, ε_minus = spectral_decomposition(ε)
    
    # 2. 计算纯净材料的拉伸能量 Ψ0+ 和压缩能量 Ψ0-
    # tr(ε) 也就是应变张量的迹
    tr_ε = tr(ε)
    
    Ψ0_plus = (mat.λ / 2) * macauley_plus(tr_ε)^2 + mat.μ * tr(ε_plus ⋅ ε_plus)
    Ψ0_minus = (mat.λ / 2) * macauley_minus(tr_ε)^2 + mat.μ * tr(ε_minus ⋅ ε_minus)
    
    # 3. 引入相场退化函数 g(d) = (1-d)^2 + k0
    # 只有拉伸能量受到损伤的削弱，压缩能量保持完整
    g_d = (1.0 - d)^2 + mat.k_tol
    
    # 4. 返回总弹性能密度
    return g_d * Ψ0_plus + Ψ0_minus
end

"""
    tensile_energy_density(ε::SymmetricTensor{2,2}, mat::PhaseFieldMaterial)

仅计算拉伸能量 Ψ0+，用于相场不可逆历史变量 H 的更新。
对应论文 Eq. (18) 中的 Ψ0+(ε(u)).
"""
function tensile_energy_density(ε::SymmetricTensor{2,2,T}, mat::PhaseFieldMaterial) where T
    ε_plus, _ = spectral_decomposition(ε)
    tr_ε = tr(ε)
    return (mat.λ / 2) * macauley_plus(tr_ε)^2 + mat.μ * tr(ε_plus ⋅ ε_plus)
end
