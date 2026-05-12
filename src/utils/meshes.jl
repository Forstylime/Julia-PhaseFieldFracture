# src/utils/meshes.jl

using Ferrite

"""
生成一个带有初始裂纹标记的 2D 正方形网格 (用于基础拉伸测试)
- L: 边长
- nx, ny: x 和 y 方向的单元数量
"""
function create_simple_tension_grid(L::Float64, nx::Int, ny::Int)
    # 1. 生成标准的四边形网格 (从 0,0 到 L,L)
    # 采用 Quadrilateral (四边形单元)，这样对于弹塑性/断裂的精度比三角形好
    grid = generate_grid(Quadrilateral, (nx, ny), Vec(0.0, 0.0), Vec(L, L))
    
    # 2. 为网格添加面集 (Face Sets)，用于施加位移边界条件
    # 底边: y = 0
    # addfacetset!(grid, "bottom", x -> x[2] ≈ 0.0)
    # 顶边: y = L
    # addfacetset!(grid, "top", x -> x[2] ≈ L)
    
    # 3. 为网格添加节点集 (Node Sets)，用于定义预制裂纹
    # 经典相场设置：我们在左侧中间水平切一刀，长度为 L/2
    # 这里我们把落在这个线段上的节点找出来，打上 "crack" 标签
    addnodeset!(grid, "crack", x -> x[2] ≈ L/2 && x[1] <= L/2)

    return grid
end