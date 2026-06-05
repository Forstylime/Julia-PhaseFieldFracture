# src/fem/assembly.jl

using Ferrite
using Tensors
using LinearAlgebra
using SparseArrays

"""
    assemble_u!(K, R, dh_u, dh_d, u, d, mat, cv_u, cv_d)

组装位移场 u 的切线刚度矩阵 K 和残差向量 R (内力)。
这是非线性牛顿迭代的核心步骤。
"""
function assemble_u!(
    K::AbstractMatrix{T}, R::AbstractVector{T}, # 允许 T 类型 (Float64 或 Dual)
    dh_u::DofHandler, dh_d::DofHandler,
    u_global::AbstractVector{T}, d_global::AbstractVector, # 位移设为 T
    mat::PhaseFieldMaterial, 
    cv_u::CellValues, cv_d::CellValues
) where T <: Real # 使用参数化类型 T
    # 初始化汇编器，自动把单元矩阵加到全局稀疏矩阵的对应位置
    assembler = start_assemble(K, R)
    
    n_basefuncs_u = getnbasefunctions(cv_u)
    
    # 获取每个单元的局部矩阵和向量缓存
    Ke = zeros(T, n_basefuncs_u, n_basefuncs_u)
    Re = zeros(T, n_basefuncs_u)
    
    # 遍历每一个网格单元 (同时遍历 u 和 d 的自由度布局)
    for (cell_u, cell_d) in zip(CellIterator(dh_u), CellIterator(dh_d))
        # 重新初始化形函数和雅可比矩阵
        reinit!(cv_u, cell_u)
        reinit!(cv_d, cell_d)
        
        # 提取当前单元上节点的全局自由度值
        u_loc = u_global[celldofs(cell_u)]
        d_loc = d_global[celldofs(cell_d)]
        
        fill!(Ke, 0.0)
        fill!(Re, 0.0)
        
        # 遍历单元内的积分点 (Gauss Quadrature points)
        for q_point in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, q_point) # 雅可比行列式乘积分权重
            
            # 计算该积分点上的应变和相场值
            ε_q = function_symmetric_gradient(cv_u, q_point, u_loc)
            d_q = function_value(cv_d, q_point, d_loc)
            
            # 能量谱分解，在特征值重复时做二阶自动微分容易产生 NaN 的问题，
            # 已在 constitutive.jl 中添加微小扰动 ε_pert 以避免这个问题。
            ψ(ε) = elastic_energy_density(ε, d_q, mat)
            
            # 直接使用 Tensors.jl 求一阶导得到应力 σ，求二阶导得到四阶材料刚度张量 ℂ
            σ = Tensors.gradient(ψ, ε_q)
            ℂ = Tensors.hessian(ψ, ε_q)
            
            # 组装单元残差和刚度矩阵
            for i in 1:n_basefuncs_u
                δε = shape_symmetric_gradient(cv_u, q_point, i)
                # 残差: 内力向量 (σ:δε)
                Re[i] += (σ ⊡ δε) * dΩ 
                
                for j in 1:n_basefuncs_u
                    Δε = shape_symmetric_gradient(cv_u, q_point, j)
                    # 切线刚度: (δε:ℂ:Δε)
                    Ke[i, j] += (δε ⊡ ℂ ⊡ Δε) * dΩ
                end
            end
        end
        # 将单元矩阵推入全局矩阵
        assemble!(assembler, celldofs(cell_u), Ke, Re)
    end
end

"""
    assemble_d!(K, F, dh_d, H, mat, cv_d)

组装相场 d 的线性方程组 K*d = F。
对应论文公式 (17) 的整理结果。
"""
function assemble_d!(
    K::SparseMatrixCSC, F::Vector{Float64}, 
    dh_d::DofHandler, H::Vector{Float64}, 
    mat::PhaseFieldMaterial, cv_d::CellValues
)
    assembler = start_assemble(K, F)
    n_basefuncs_d = getnbasefunctions(cv_d)
    
    Ke = zeros(n_basefuncs_d, n_basefuncs_d)
    Fe = zeros(n_basefuncs_d)
    
    qp_count = 1
    for cell in CellIterator(dh_d)
        reinit!(cv_d, cell)
        
        fill!(Ke, 0.0)
        fill!(Fe, 0.0)
        
        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)
            H_q = H[qp_count] # 取出该积分点最新的历史变量
            
            # 将弱形式 (17) 重排为 A*d = B 的形式
            coef_d    = mat.gc / mat.l + 2.0 * H_q
            coef_grad = mat.gc * mat.l
            
            for i in 1:n_basefuncs_d
                δd  = shape_value(cv_d, q_point, i)
                ∇δd = shape_gradient(cv_d, q_point, i)
                
                # 右端项载荷向量
                Fe[i] += (2.0 * H_q * δd) * dΩ
                
                for j in 1:n_basefuncs_d
                    Δd  = shape_value(cv_d, q_point, j)
                    ∇Δd = shape_gradient(cv_d, q_point, j)
                    
                    # 刚度矩阵
                    Ke[i, j] += (coef_d * δd * Δd + coef_grad * (∇δd ⋅ ∇Δd)) * dΩ
                end
            end
            qp_count += 1
        end
        assemble!(assembler, celldofs(cell), Ke, Fe)
    end
end

"""
    assemble_mass_matrix_d!(M, dh_d, cv_d)

组装相场 d 的全局一致质量矩阵 M_ij = ∫_Ω N_i·N_j dΩ。
该矩阵与材料参数和状态变量无关，网格拓扑不变时为常数矩阵，
只需在求解器初始化时组装一次。
"""
function assemble_mass_matrix_d!(
    M::SparseMatrixCSC, dh_d::DofHandler, cv_d::CellValues,
)
    assembler = start_assemble(M)
    n_basefuncs_d = getnbasefunctions(cv_d)
    Me = zeros(n_basefuncs_d, n_basefuncs_d)

    for cell in CellIterator(dh_d)
        reinit!(cv_d, cell)
        fill!(Me, 0.0)

        for q_point in 1:getnquadpoints(cv_d)
            dΩ = getdetJdV(cv_d, q_point)

            for i in 1:n_basefuncs_d
                N_i = shape_value(cv_d, q_point, i)
                for j in 1:n_basefuncs_d
                    N_j = shape_value(cv_d, q_point, j)
                    Me[i, j] += N_i * N_j * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Me)
    end
end

"""
使用统一的 DofHandler 组装整体刚度矩阵 K_aa 和残差 r_a。
这就是论文 Eq. 23 和 Eq. 24 中的 K_aa 和 f_int。
"""
function assemble_monolithic!(
    K_aa::Union{AbstractMatrix{T}, Nothing}, r_a::AbstractVector{T},
    dh_a::DofHandler, x_global::AbstractVector{T},
    H_old::AbstractVector,
    mat::PhaseFieldMaterial,
    cv_u::CellValues, cv_d::CellValues
) where T <: Real

    # 仅在需要装配刚度矩阵时初始化装配器
    assembler = K_aa !== nothing ? start_assemble(K_aa, r_a) : nothing

    # 获取 u 和 d 在单元矩阵中的索引范围
    u_range = Ferrite.dof_range(dh_a, :u)
    d_range = Ferrite.dof_range(dh_a, :d)

    n_dofs = ndofs_per_cell(dh_a)
    Ke = zeros(T, n_dofs, n_dofs)
    Re = zeros(T, n_dofs)

    qp_count = 1

    for cell in CellIterator(dh_a)
        reinit!(cv_u, cell)
        reinit!(cv_d, cell)

        c_dofs = celldofs(cell)
        x_loc = x_global[c_dofs]
        u_loc = x_loc[u_range]
        d_loc = x_loc[d_range]

        fill!(Re, zero(T))
        if K_aa !== nothing
            fill!(Ke, zero(T))
        end

        for q_point in 1:getnquadpoints(cv_u)
            dΩ = getdetJdV(cv_u, q_point)

            ε_q  = function_symmetric_gradient(cv_u, q_point, u_loc)

            d_q  = function_value(cv_d, q_point, d_loc)
            ∇d_q = function_gradient(cv_d, q_point, d_loc)

            # ---- 退化函数 ----
            g_q  = (1.0 - d_q)^2 + mat.k_tol
            dg_q = -2.0 * (1.0 - d_q)

            # ---- 分解应力 ----
            split_results = miehe_spectral_decomposition(ε_q, T(mat.λ), T(mat.μ))
            σ_real = g_q * split_results.σ

            # ---- 历史变量 ----
            H_old_q = H_old[qp_count]
            is_active = split_results.ψ_pos > H_old_q
            H_q = is_active ? split_results.ψ_pos : H_old_q

            # ---- 相场系数 ----
            coef_d    = mat.gc / mat.l + 2.0 * H_q
            coef_grad = mat.gc * mat.l

            # ----------------------------------------------------
            # 填入单元残差 Re（及可选的 Ke）
            # ----------------------------------------------------
            # 1. 位移场部分
            for i in eachindex(u_range)
                I_u = u_range[i]
                δε_i = shape_symmetric_gradient(cv_u, q_point, i)

                # R_u: 内力残差
                Re[I_u] += (σ_real ⊡ δε_i) * dΩ

                if K_aa !== nothing
                    ℂ_damage = g_q * split_results.ℂ
                    # K_uu: ∂R_u/∂u
                    for j in eachindex(u_range)
                        J_u = u_range[j]
                        Δε_j = shape_symmetric_gradient(cv_u, q_point, j)
                        Ke[I_u, J_u] += (δε_i ⊡ ℂ_damage ⊡ Δε_j) * dΩ
                    end

                    # K_ud: ∂R_u/∂d
                    for j in eachindex(d_range)
                        J_d = d_range[j]
                        N_j = shape_value(cv_d, q_point, j)
                        Ke[I_u, J_d] += dg_q * (split_results.σ ⊡ δε_i) * N_j * dΩ
                    end
                end
            end

            # 2. 相场部分
            for i in eachindex(d_range)
                I_d = d_range[i]
                δd_i  = shape_value(cv_d, q_point, i)
                ∇δd_i = shape_gradient(cv_d, q_point, i)

                # R_d: 相场残差
                Re[I_d] += coef_d * δd_i * d_q * dΩ
                Re[I_d] += coef_grad * (∇δd_i ⋅ ∇d_q) * dΩ
                Re[I_d] -= 2.0 * H_q * δd_i * dΩ

                if K_aa !== nothing
                    # K_du: ∂R_d/∂u
                    if is_active
                        for j in eachindex(u_range)
                            J_u = u_range[j]
                            Δε_j = shape_symmetric_gradient(cv_u, q_point, j)
                            Ke[I_d, J_u] += 2.0 * (d_q - 1.0) * δd_i * (split_results.σ ⊡ Δε_j) * dΩ
                        end
                    end

                    # K_dd: ∂R_d/∂d
                    for j in eachindex(d_range)
                        J_d = d_range[j]
                        Δd_j  = shape_value(cv_d, q_point, j)
                        ∇Δd_j = shape_gradient(cv_d, q_point, j)
                        Ke[I_d, J_d] += (coef_d * δd_i * Δd_j + coef_grad * (∇δd_i ⋅ ∇Δd_j)) * dΩ
                    end
                end
            end

            qp_count += 1
        end

        # 区分装配模式
        if K_aa !== nothing
            assemble!(assembler, c_dofs, Ke, Re)
        else
            r_a[c_dofs] .+= Re  # 绕过装配器，直接原位累加（支持自动微分）
        end
    end
end