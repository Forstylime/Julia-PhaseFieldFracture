# src/physics/constitutive.jl

using Tensors
using LinearAlgebra

"""
    PhaseFieldMaterial

存储线弹性相场断裂模型的材料参数。
论文采用 Lamé 常数 λ 和 μ，以及断裂能 gc 和尺度参数 l。
"""
Base.@kwdef struct PhaseFieldMaterial
    E::Float64        = 25840   # 杨氏模量 (MPa)
    ν::Float64        = 0.18    # 泊松比
    # Lamé 常数 (根据论文，采用平面应力，而不是平面应变公式)
    λ::Float64        = (E * ν) / (1 - ν^2) # (E * ν) / ((1 + ν) * (1 - 2ν)) #
    μ::Float64        = E / (2 * (1 + ν))
    # 相场断裂参数
    gc::Float64       = 0.65    # 临界断裂能 (N/mm)
    l::Float64        = 10      # 长度尺度参数 (mm)
    k_tol::Float64    = 1e-8    # 残余刚度参数 k0 (防止 d=1 时刚度矩阵奇异)
end

"""
    macauley_plus(x) / macauley_minus(x)

Macaulay 括号，用于提取正数或负数部分，对应论文公式里的 ⟨x⟩±。
"""
@inline macauley_plus(x::Real) = max(x, zero(x))
@inline macauley_minus(x::Real) = min(x, zero(x))

"""
    strain_spectral_split(ε::SymmetricTensor{2, 2})

对 2D 应变张量进行谱分解，返回正部分 ε+ 和负部分 ε-。
对应论文 Eq. (4)。
"""
function strain_spectral_split(ε::SymmetricTensor{2, 2, T}) where T
    # 对对称张量求特征值和特征向量
    # 注意：因为 ForwardDiff 传进来的类型 T 可能是 Dual 数，
    # 这里我们保证类型泛型 T，以支持自动微分穿透！

    # 如果应变几乎为 0，直接返回零张量，避免特征值重根引发的 AD 求导 NaN
    if norm(ε) < 1e-14
        return zero(ε), zero(ε)
    end
    
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
    ε_pert = SymmetricTensor{2, 2}((1e-14, 0.0, -1e-14))
    ε = ε + ε_pert # 添加微小扰动，避免特征值重根引发的 AD 求导 NaN
    ε_plus, ε_minus = strain_spectral_split(ε)
    
    # 2. 计算纯净材料的拉伸能量 Ψ0+ 和压缩能量 Ψ0-
    # tr(ε) 也就是应变张量的迹
    tr_ε = tr(ε)
    
    Ψ0_plus = (mat.λ / 2) * macauley_plus(tr_ε)^2 + mat.μ * tr(ε_plus ⋅ ε_plus) # 也可以写 mat.μ * (ε_plus ⊡ ε_plus)
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
    ε_pert = SymmetricTensor{2, 2}((1e-14, 0.0, -1e-14))
    ε = ε + ε_pert # 添加微小扰动，避免特征值重根引发的 AD 求导 NaN
    ε_plus, _ = strain_spectral_split(ε)
    tr_ε = tr(ε)
    return (mat.λ / 2) * macauley_plus(tr_ε)^2 + mat.μ * tr(ε_plus ⋅ ε_plus)
end



"""
根据 Miehe 的谱分解，计算应力、拉伸能和一致切线刚度张量。
正确处理了主应变相等的退化情况。

返回: NamedTuple (σ=..., ψ_pos=..., ℂ=...)
"""
function miehe_spectral_decomposition(
    ε::SymmetricTensor{2, 2, T}, 
    λ_bar::T, 
    μ::T
) where T
    
    I2 = one(SymmetricTensor{2, 2, T})
    I4_sym = one(SymmetricTensor{4, 2, T}) # 4阶对称单位张量

    # --- 1. 特征值分解 ---
    vals, vecs = eigen(ε)
    λ₁, λ₂ = vals[1], vals[2]
    #n₁, n₂ = vecs[1], vecs[2]
    n₁ = vecs[:, 1]
    n₂ = vecs[:, 2]

    # --- 2. 计算拉伸相关量 ---
    pos(x) = x > 0.0 ? x : 0.0
    H(x) = x > 0.0 ? 1.0 : 0.0

    λ₁⁺, λ₂⁺ = pos(λ₁), pos(λ₂)
    
    tr_ε = λ₁ + λ₂
    tr_ε_pos = pos(tr_ε)

    ε_pos = λ₁⁺ * (n₁ ⊗ n₁) + λ₂⁺ * (n₂ ⊗ n₂)
    ψ_pos = (λ_bar / 2.0) * tr_ε_pos^2 + μ * (λ₁⁺^2 + λ₂⁺^2)
    σ = λ_bar * tr_ε_pos * I2 + 2.0 * μ * ε_pos

    # --- 3. 计算一致切线刚度张量 ℂ ---
    ℂ = zero(SymmetricTensor{4, 2, T})

    # 3.1 体积部分 (来自 tr(ε))
    ℂ += λ_bar * H(tr_ε) * (I2 ⊗ I2)

    # 3.2 偏量部分的主方向贡献
    #P₁ = n₁ ⊗ n₁
    #P₂ = n₂ ⊗ n₂
    P₁ = symmetric(n₁ ⊗ n₁)
    P₂ = symmetric(n₂ ⊗ n₂)
    ℂ += 2.0 * μ * H(λ₁) * (P₁ ⊗ P₁)
    ℂ += 2.0 * μ * H(λ₂) * (P₂ ⊗ P₂)

    # 3.3 偏量部分的耦合项（最复杂的部分）
    if !isapprox(λ₁, λ₂; atol=1e-8, rtol=1e-6)
        # 非退化情况
        coef = 2.0 * μ * (λ₁⁺ - λ₂⁺) / (λ₁ - λ₂)
        # 构造耦合张量 M₁₂
        M₁₂ = (o_outer(n₁, n₂, n₂, n₁) + o_outer(n₁, n₂, n₁, n₂) +
               o_outer(n₂, n₁, n₂, n₁) + o_outer(n₂, n₁, n₁, n₂)) / 2.0
        ℂ += coef * M₁₂
    else
        # 退化情况 (λ₁ ≈ λ₂), 此时分子分母的导数之比为 H(λ₁)
        coef = 2.0 * μ * H(λ₁)
        # I₄_dev = I₄_sym - (1/2)*(I2⊗I2), P₁ + P₂ = I₂
        # 最终退化为标准的各向同性偏量切线
        ℂ += coef * (I4_sym - P₁ ⊗ P₁ - P₂ ⊗ P₂) 
    end
    
    return (
        σ_pos = σ,        # 也就是之前的 σ
        ℂ_pos = ℂ,        # 也就是之前的 ℂ
        ψ_pos = ψ_pos
    )
end

# 辅助函数 o_outer(a,b,c,d) = a⊗b⊗c⊗d
function o_outer(a, b, c, d)
    symmetric((a ⊗ b) ⊗ (c ⊗ d))
end