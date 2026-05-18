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
    K::SparseMatrixCSC, R::Vector{Float64}, 
    dh_u::DofHandler, dh_d::DofHandler,
    u_global::Vector{Float64}, d_global::Vector{Float64}, 
    mat::PhaseFieldMaterial, 
    cv_u::CellValues, cv_d::CellValues
)
    # 初始化汇编器，自动把单元矩阵加到全局稀疏矩阵的对应位置
    assembler = start_assemble(K, R)
    
    n_basefuncs_u = getnbasefunctions(cv_u)
    
    # 获取每个单元的局部矩阵和向量缓存
    Ke = zeros(n_basefuncs_u, n_basefuncs_u)
    Re = zeros(n_basefuncs_u)
    
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
            
            # 谱分解能量在特征值重复时做二阶自动微分容易产生 NaN。
            # 位移方程使用稳定的退化线弹性切线，历史变量仍在 `update_history!`
            # 中通过拉伸谱能量控制裂纹不可逆演化。
            g_q = (1.0 - d_q)^2 + mat.k_tol
            ψ(ε) = g_q * ((mat.λ / 2) * tr(ε)^2 + mat.μ * tr(ε ⋅ ε))
            
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
