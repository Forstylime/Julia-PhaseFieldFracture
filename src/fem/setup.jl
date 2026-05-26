"""
    SquareTensionSetup

保存方形拉伸相场断裂算例的有限元初始化结果。也可以用于其他几何相似的算例（如 L 形拉伸），只要保证网格和边界条件设置与求解器要求一致。

该结构体把后续求解器需要反复使用的对象集中在一起：计算网格、位移场和
相场的自由度处理器、两类场的约束处理器，以及用于设置预制裂纹的节点编号。
"""
Base.@kwdef struct SquareTensionSetup{G,DHU,DHD,CHU,CHD,N}
    grid::G
    dh_u::DHU
    dh_d::DHD
    ch_u::CHU
    ch_d::CHD
    crack_nodes::N
    final_displacement::Float64
end

"""
    make_square_tension_grid(cells = (100, 100); xlims = (-1.0, 1.0), ylims = (-1.0, 1.0),
                             refine_y_center = nothing, refine_ratio = 3.0)

生成方形拉伸算例使用的二维四边形结构网格，支持局部细化。

# 参数
- `cells`: x、y 方向的单元数量。
- `xlims`, `ylims`: 计算域的坐标范围。
- `refine_y_center`: y 方向加密中心。`nothing` 时 y 方向均匀划分。
  设为 `0.0` 可在裂纹高度处加密 y 方向网格。
- `refine_ratio`: 加密强度，> 1 时中心更密（建议 3.0 ~ 5.0），
  = 1 时退化为均匀分布。

# 使用方法
```julia
# 均匀网格（原有行为）
grid = make_square_tension_grid((100, 100))

# y 方向在 y=0（裂纹线）附近加密
grid = make_square_tension_grid((100, 100); refine_y_center = 0.0, refine_ratio = 4.0)
```
"""
function make_square_tension_grid(
    cells::NTuple{2,Int} = (100, 100);
    xlims::NTuple{2,Float64} = (-1.0, 1.0),
    ylims::NTuple{2,Float64} = (-1.0, 1.0),
    refine_y_center::Union{Float64,Nothing} = 0.0,
    refine_x_center::Union{Float64,Nothing} = nothing,
    refine_ratio::Float64 = 3.0,
)
    nx, ny = cells
    n_nodes_x, n_nodes_y = nx + 1, ny + 1

    # 先用 Ferrite 生成均匀网格，获得正确的拓扑和边界 facet set
    grid = Ferrite.generate_grid(
        Ferrite.Quadrilateral, (nx, ny),
        Ferrite.Vec(xlims[1], ylims[1]),
        Ferrite.Vec(xlims[2], ylims[2]),
    )

    # 如果不需要加密，直接返回均匀网格
    if refine_x_center === nothing && refine_y_center === nothing
        return grid
    end

    # x 方向坐标（n_nodes_x 个节点）
    if refine_x_center !== nothing
        Lx = xlims[2] - xlims[1]
        xc = refine_x_center
        x = refine_range(Lx, n_nodes_x; center = xc, ratio = refine_ratio)
        x .+= (xc - (x[1] + x[end]) / 2)
        xv = collect(x)
    else
        xv = collect(range(Float64(xlims[1]), Float64(xlims[2]), length = n_nodes_x))
    end

    # y 方向坐标（n_nodes_y 个节点）
    if refine_y_center !== nothing
        Ly = ylims[2] - ylims[1]
        yc = refine_y_center
        y = refine_range(Ly, n_nodes_y; center = yc, ratio = refine_ratio)
        y .+= (yc - (y[1] + y[end]) / 2)
        yv = collect(y)
    else
        yv = collect(range(Float64(ylims[1]), Float64(ylims[2]), length = n_nodes_y))
    end

    # 修改节点坐标以匹配非均匀分布
    refine_grid!(grid, xv, yv)
    return grid
end

"""
    create_l_shape_grid(msh_file = "data/mesh/l_shape.msh")
    生成 L 形算例使用的二维四边形结构网格。
    网格一般已在Gmsh中生成好，直接从 .msh 文件读取。
"""
function create_l_shape_grid(msh_file = "data/mesh/l_shape.msh")
    cache_file = msh_file * ".jls"
    if isfile(cache_file) && mtime(cache_file) >= mtime(msh_file)
        println("从缓存加载网格: ", cache_file)
        return deserialize(cache_file)
    end
    println("解析 .msh 文件: ", msh_file)
    grid = FerriteGmsh.togrid(msh_file)
    serialize(cache_file, grid)
    println("网格已缓存至: ", cache_file)
    return grid
end

"""
    create_staggered_dofhandlers(grid)

为交错求解格式创建位移场和相场各自的自由度处理器。

# 功能
- 位移场 `:u` 使用一阶 Lagrange 四边形插值，并通过 `^2` 表示二维向量场。
- 相场 `:d` 使用一阶 Lagrange 四边形标量插值。
- 两个场使用独立的 `DofHandler`，便于交错求解器分别组装和求解。

# 使用方法
```julia
dh_u, dh_d = create_staggered_dofhandlers(grid)
```
"""
function create_staggered_dofhandlers(grid)
    # 位移自由度处理器：每个节点有两个分量，对应 ux 和 uy。
    dh_u = Ferrite.DofHandler(grid)
    ip_u = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2
    Ferrite.add!(dh_u, :u, ip_u)
    # `close!` 会冻结自由度布局并建立单元到全局自由度的映射。
    Ferrite.close!(dh_u)

    # 相场自由度处理器：每个节点一个标量损伤变量 d。
    dh_d = Ferrite.DofHandler(grid)
    ip_d = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()
    Ferrite.add!(dh_d, :d, ip_d)
    # 关闭后才能查询自由度数量、组装矩阵或创建约束。
    Ferrite.close!(dh_d)

    println("位移自由度数量: ", Ferrite.ndofs(dh_u))
    println("相场自由度数量: ", Ferrite.ndofs(dh_d))

    return dh_u, dh_d
end

"""
    create_displacement_constraints(dh_u, grid; final_displacement = 0.0)

创建位移场的 Dirichlet 边界条件。

# 功能
- 底边 `"bottom"`或者顶边 `"top"` 的两个位移分量均固定为 0，消除刚体运动。
- 顶边 `"top"`或者右边界 `"right"` 的竖向位移分量施加为 `t * final_displacement`，
  其中 `t` 由 `Ferrite.update!(ch_u, t)` 控制，可用于增量加载。
- 水平位移分量在顶边不额外约束。

# 使用方法
```julia
ch_u = create_displacement_constraints(dh_u, grid; final_displacement = 0.01)
```
"""
function create_displacement_constraints(dh_u, grid; final_displacement = 0.0)
    ch_u = Ferrite.ConstraintHandler(dh_u)

    # 读取网格生成阶段创建的边界 facet set。
    top = Ferrite.getfacetset(grid, "top")
    right = Ferrite.getfacetset(grid, "right")

    # 顶边 ux、uy 全固定，作为拉伸试样的支承边界。
    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, top, (x, t) -> zeros(2), [1, 2]),
    )
    # 右边界只约束第 2 个位移分量 uy；加载幅值通过时间/载荷参数 t 缩放。
    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, right, (x, t) -> t * final_displacement, 2),
    )

    # 关闭约束处理器并在 t = 0 时初始化约束值。
    Ferrite.close!(ch_u)
    Ferrite.update!(ch_u, 0.0)
    return ch_u
end

"""
    create_phase_field_constraints(dh_d, crack_nodes)

创建相场变量的约束处理器。

# 功能
- 在 `crack_nodes` 上施加 Dirichlet BC，固定 d = 1 (预制裂纹)。
- 若 `crack_nodes` 为空，则不添加任何约束（相场完全自由演化）。
- 返回已关闭并初始化过的 `ConstraintHandler`。

# 使用方法
```julia
ch_d = create_phase_field_constraints(dh_d, crack_nodes)
```
"""
function create_phase_field_constraints(dh_d, crack_nodes)
    ch_d = Ferrite.ConstraintHandler(dh_d)
    if !isempty(crack_nodes)
        Ferrite.add!(ch_d, Ferrite.Dirichlet(:d, collect(crack_nodes), (x, t) -> 1.0))
    end
    Ferrite.close!(ch_d)
    Ferrite.update!(ch_d, 0.0)
    return ch_d
end

"""
    initial_crack_nodes(grid; y = 0.0, x_min = -1.0, x_max = 0.0, half_width = 1e-10)

查找位于预制裂纹线段附近的网格节点。

# 功能
- 在节点坐标中筛选满足 `x_min <= x <= x_max` 的节点。
- 同时要求节点的 `y` 坐标与给定裂纹高度 `y` 的距离不超过 `half_width`。
- 返回节点编号向量，供求解器或初始化过程把这些节点设置为已损伤/裂纹状态。

# 使用方法
```julia
crack_nodes = initial_crack_nodes(grid; y = 0.0, x_min = -1.0, x_max = 0.0)
```
"""
function initial_crack_nodes(
    grid;
    y = 0.0,
    x_min = -1.0,
    x_max = 0.0,
    half_width = 1e-8,
)
    nodes = Int[]
    # 遍历 Ferrite 网格节点；`pairs` 同时给出节点编号和节点对象。
    for (i, node) in pairs(grid.nodes)
        xcoord = node.x[1]
        ycoord = node.x[2]
        # 用一个很小的厚度容差捕捉裂纹线，避免浮点坐标比较过于脆弱。
        if x_min <= xcoord <= x_max && abs(ycoord - y) <= half_width
            push!(nodes, i)
        end
    end
    return nodes
end

"""
    setup_square_tension(; cells = (50, 50), final_displacement = 0.0,
                           crack_y = 0.0, crack_x_min = -1.0,
                           crack_x_max = 0.0, crack_half_width = 1e-10)

一站式构建方形拉伸相场断裂算例的有限元初始化对象。

# 脚本功能
本文件提供方形拉伸问题的前处理工具：生成网格、建立位移场/相场自由度、
设置拉伸边界条件、定位预制裂纹节点，并把这些对象封装到
`SquareTensionSetup` 中。求解脚本可以直接使用返回对象进入装配和求解阶段。

# 使用方法
```julia
using PhaseFieldFracture

setup = setup_square_tension(
    cells = (50, 50),
    final_displacement = 0.01,
    crack_y = 0.0,
    crack_x_min = -1.0,
    crack_x_max = 0.0,
)

grid = setup.grid
dh_u = setup.dh_u
dh_d = setup.dh_d
crack_nodes = setup.crack_nodes
```

# 参数说明
- `cells`：x、y 方向的单元数量。
- `final_displacement`：顶边竖向位移加载幅值。
- `crack_y`：预制裂纹所在的 y 坐标。
- `crack_x_min` / `crack_x_max`：预制裂纹在线段方向上的 x 坐标范围。
- `crack_half_width`：筛选裂纹节点时使用的半宽容差。
"""
function setup_square_tension(;
    cells::NTuple{2,Int} = (100, 100),
    final_displacement = 0.0,
    crack_y = 0.0,
    crack_x_min = -1.0,
    crack_x_max = 0.0,
    crack_half_width = 1e-10,
)
    # 1. 创建计算网格，后续所有自由度和边界集合都基于同一个 grid。
    grid = make_square_tension_grid(cells)
    # 2. 分别为位移场 u 和相场 d 建立自由度编号，适配交错求解流程。
    dh_u, dh_d = create_staggered_dofhandlers(grid)
    # 3. 设置力学边界条件：底边固定，顶边施加竖向位移加载。
    ch_u = create_displacement_constraints(dh_u, grid; final_displacement)
    # 4. 根据几何坐标筛选预制裂纹节点，供初始损伤场使用。
    crack_nodes = initial_crack_nodes(
        grid;
        y = crack_y,
        x_min = crack_x_min,
        x_max = crack_x_max,
        half_width = crack_half_width,
    )
    # 5. 创建相场约束处理器，在裂纹节点上固定 d = 1。
    ch_d = create_phase_field_constraints(dh_d, crack_nodes)

    # 统一封装初始化结果，减少求解脚本需要手动传递的对象数量。
    return SquareTensionSetup(; grid, dh_u, dh_d, ch_u, ch_d, crack_nodes, final_displacement)
end

"""
    setup_l_tension(; msh_file = "data/mesh/l_shape.msh", final_displacement = 0.0)
    构建 L 形拉伸相场断裂算例的有限元初始化对象。
    一般不需要预置裂纹，L 形几何本身会引导裂纹从内角处自然萌生。
"""
function setup_l_tension(;
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = 0.0,
)
    # 1. 创建计算网格，后续所有自由度和边界集合都基于同一个 grid。
    grid = create_l_shape_grid(msh_file)
    # 2. 分别为位移场 u 和相场 d 建立自由度编号，适配交错求解流程。
    dh_u, dh_d = create_staggered_dofhandlers(grid)
    # 3. 设置力学边界条件：上边(top)固定，右边(right)施加竖直向下位移加载。
    ch_u = create_displacement_constraints(dh_u, grid; final_displacement)
    # 4. 创建相场约束处理器。
    ch_d = create_phase_field_constraints(dh_d, Int[]) # L 形算例通常不需要预置裂纹，传入空节点列表。

    # 统一封装初始化结果，减少求解脚本需要手动传递的对象数量。
    return SquareTensionSetup(; grid, dh_u, dh_d, ch_u, ch_d, crack_nodes = Int[], final_displacement)
end