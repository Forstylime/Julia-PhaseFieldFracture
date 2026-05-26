"""
    计算各类能量的函数
"""

"""
    surface energy
    计算单元的表面能，积分点上的表面能密度为 (g_c / (2*l)) * (d^2 + l^2 * |∇d|^2)。
    其中 g_c 是材料的断裂韧性，l 是长度尺度参数，d 是相场变量，∇d 是相场的梯度。
    参考论文中的 Eq. (5)
"""
function surface_energy(mat::PhaseFieldMaterial, cv_d::CellValues)
    ∇d = shape_gradient(cv_d)
    dΩ = getdetJdV(cv_d)

    # 计算积分点上的表面能密度 (g_c / (2*l)) * (d^2 + l^2 * |∇d|^2)
    energy_density = (mat.gc / (2 * mat.l)) * (shape_value(cv_d, :)^2 + mat.l^2 * sum(∇d.^2, dims=1))
    
    # 对所有积分点求和得到单元的表面能
    return sum(energy_density) * dΩ
end