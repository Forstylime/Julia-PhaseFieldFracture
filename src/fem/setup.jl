"""
    TensionSetup

保存方形拉伸相场断裂算例的有限元初始化结果。也可以用于其他几何相似的算例（如 L 形拉伸），只要保证网格和边界条件设置与求解器要求一致。

该结构体把后续求解器需要反复使用的对象集中在一起：计算网格、位移场和
相场的自由度处理器、两类场的约束处理器，以及用于设置预制裂纹的节点编号。
"""
Base.@kwdef struct TensionSetup{G,DHU,DHD,CHU,CHD,N}
    grid::G
    dh_u::DHU
    dh_d::DHD
    ch_u::CHU
    ch_d::CHD
    crack_nodes::N
    final_displacement::Float64
end

"""
    create_grid(msh_file = "data/mesh/l_shape.msh")
    生成 L 形算例使用的二维四边形结构网格。
    网格一般已在Gmsh中生成好，直接从 .msh 文件读取。
"""
function create_grid(msh_file = "data/mesh/l_shape.msh")
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
function create_displacement_constraints_square(dh_u, grid; final_displacement = 0.0)
    ch_u = Ferrite.ConstraintHandler(dh_u)

    # 读取网格生成阶段创建的边界 facet set。
    left = Ferrite.getfacetset(grid, "left")
    right = Ferrite.getfacetset(grid, "right")

    # 左边界 ux、uy 全固定，作为拉伸试样的支承边界。
    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, left, (x, t) -> zeros(2), [1, 2]),
    )
    # 右边界只约束第 1 个位移分量 ux；加载幅值通过时间/载荷参数 t 缩放。
    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, right, (x, t) -> t * final_displacement, 1),
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
    initial_crack_nodes(grid; ...)

查找位于预制裂纹线段附近的网格节点。

# 功能
- 在节点坐标中筛选满足条件的节点。
- 同时要求节点的坐标与给定裂纹位置的距离不超过 `half_width`。
- 返回节点编号向量，供求解器或初始化过程把这些节点设置为已损伤/裂纹状态。
```
"""
function initial_crack_nodes(
    grid;
    x = 1.0,
    y_min = 0.0,
    y_max = 1.0,
    half_width = 1e-8,
)
    nodes = Int[]
    # 遍历 Ferrite 网格节点；`pairs` 同时给出节点编号和节点对象。
    for (i, node) in pairs(grid.nodes)
        xcoord = node.x[1]
        ycoord = node.x[2]
        # 用一个很小的厚度容差捕捉裂纹线，避免浮点坐标比较过于脆弱。
        if y_min <= ycoord <= y_max && abs(xcoord - x) <= half_width
            push!(nodes, i)
        end
    end
    return nodes
end

"""
    setup_square_tension(; msh_file = "data/mesh/square.msh", final_displacement = 0.0)

一站式构建方形拉伸相场断裂算例的有限元初始化对象。

# 脚本功能
本文件提供方形拉伸问题的前处理工具：生成网格、建立位移场/相场自由度、
设置拉伸边界条件、定位预制裂纹节点，并把这些对象封装到
`TensionSetup` 中。求解脚本可以直接使用返回对象进入装配和求解阶段。

# 参数说明
- `msh_file`：Gmsh 网格文件路径。
- `final_displacement`：顶边竖向位移加载幅值。
"""
function setup_square_tension(;
    msh_file = "data/mesh/square.msh",
    final_displacement = 0.0,
)
    # 1. 创建计算网格，后续所有自由度和边界集合都基于同一个 grid。
    grid = create_square_tension_grid(msh_file)
    # 2. 分别为位移场 u 和相场 d 建立自由度编号，适配交错求解流程。
    dh_u, dh_d = create_staggered_dofhandlers(grid)
    # 3. 设置力学边界条件：底边固定，顶边施加竖向位移加载。
    ch_u = create_displacement_constraints_square(dh_u, grid; final_displacement)
    # 5. 创建相场约束处理器
    ch_d = create_phase_field_constraints(dh_d, Int[]) # CT 试件通过预制缺口实现裂纹自然产生，也不需要预制裂纹

    # 统一封装初始化结果，减少求解脚本需要手动传递的对象数量。
    return TensionSetup(; grid, dh_u, dh_d, ch_u, ch_d, crack_nodes = Int[], final_displacement)
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
    grid = create_grid(msh_file)
    # 2. 分别为位移场 u 和相场 d 建立自由度编号，适配交错求解流程。
    dh_u, dh_d = create_staggered_dofhandlers(grid)
    # 3. 设置力学边界条件：上边(top)固定，右边(right)施加竖直向下位移加载。
    ch_u = create_displacement_constraints(dh_u, grid; final_displacement)
    # 4. 创建相场约束处理器。
    ch_d = create_phase_field_constraints(dh_d, Int[]) # L 形算例通常不需要预置裂纹，传入空节点列表。

    # 统一封装初始化结果，减少求解脚本需要手动传递的对象数量。
    return TensionSetup(; grid, dh_u, dh_d, ch_u, ch_d, crack_nodes = Int[], final_displacement)
end



"""
|===========================================================|
|    MonolithicTensionSetup                                 |
|                                                           |
|保存算例的有限元初始化结果（整体式版本）。                 |
|只保留一个统一的 DofHandler 和 ConstraintHandler。         |
|===========================================================|
"""
struct MonolithicTensionSetup
    dir::Int
    grid::Grid
    dh::DofHandler
    ch_ref::ConstraintHandler   # 专门用于预测步
    ch_zero::ConstraintHandler  # 专门用于校正步
    final_displacement::Float64 
end

"""
    create_monolithic_dofhandler(grid)

为整体式求解格式创建一个包含位移场(:u)和相场(:d)的统一自由度处理器。
"""
function create_monolithic_dofhandler(grid)
    dh_a = Ferrite.DofHandler(grid)
    
    # 位移场：2维向量
    ip_u = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2
    Ferrite.add!(dh_a, :u, ip_u)
    
    # 相场：标量
    ip_d = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()
    Ferrite.add!(dh_a, :d, ip_d)
    
    Ferrite.close!(dh_a)

    println("总自由度数量 (u + d): ", Ferrite.ndofs(dh_a))
    return dh_a
end

"""
    create_arc_length_bcs(dh, grid, crack_nodes)

为整体式 DofHandler 创建边界条件。
"""
function create_arc_length_bcs(dh, grid, fixed_face = "top", final_displacement = -0.8, dir = 2)
    fixed = Ferrite.getfacetset(grid, fixed_face)
    right = Ferrite.getfacetset(grid, "right")

    # ==============================================================
    # 1. ch_ref: 用于预测步 (求解参考位移 Δa_λ / u_T)
    # ==============================================================
    ch_ref = Ferrite.ConstraintHandler(dh)
    
    # 固定 (0.0)
    Ferrite.add!(ch_ref, Ferrite.Dirichlet(:u, fixed, (x, t) -> zeros(2), [1, 2]))
    
    # 【Trick 核心】右侧施加竖向参考位移 1.0 (之后由 solver 里的 λ 自动放大)
    Ferrite.add!(ch_ref, Ferrite.Dirichlet(:u, right, (x, t) -> final_displacement, dir))
    
    Ferrite.close!(ch_ref)
    Ferrite.update!(ch_ref, 0.0) # 时间参数用不上了，随便传个0

    # ==============================================================
    # 2. ch_zero: 用于校正步 (求解残余力对应的位移修正 Δa_r)
    # ==============================================================
    ch_zero = Ferrite.ConstraintHandler(dh)
    
    # 固定 (0.0)
    Ferrite.add!(ch_zero, Ferrite.Dirichlet(:u, fixed, (x, t) -> zeros(2), [1, 2]))
    
    # 右侧受控端：必须锁定为 0.0 (残余修正不能改变已经由 λ 决定的位移)
    Ferrite.add!(ch_zero, Ferrite.Dirichlet(:u, right, (x, t) -> 0.0, dir))
    
    Ferrite.close!(ch_zero)
    Ferrite.update!(ch_zero, 0.0)

    return ch_ref, ch_zero
end

"""
    setup_tension_monolithic(...)

一站式初始化函数（整体式版本示例）。
"""
function setup_tension_monolithic(;
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8, # 最终的位移量
    dir = 2, # 位移方向, 1 -> x , 2 -> y
    fixed_face = "top", # 固定边
)
    grid = create_grid(msh_file)
    
    # 1. 创建整体式 DofHandler
    dh = create_monolithic_dofhandler(grid)
    
    # 2. 获取弧长法专用的两套 ConstraintHandler
    ch_ref, ch_zero = create_arc_length_bcs(dh, grid, fixed_face, final_displacement, dir)
    
    return MonolithicTensionSetup(dir, grid, dh, ch_ref, ch_zero, final_displacement)
end

"""
    MonolithicTensionSetup_MEM (适配 M-EM 算法)

保存算例的有限元初始化结果。
包含了全局单片 dh，纯相场 dh_d，以及时间依赖的 ch_a 和齐次的 ch_zero。
"""
struct MonolithicTensionSetup_MEM
    grid::Grid
    dh::DofHandler
    dh_d::DofHandler            # 【新增】：专用于组装相场质量矩阵 M_d
    ch_a::ConstraintHandler     # 【修改】：随时间 t 动态更新的边界条件
    ch_zero::ConstraintHandler  # 用于 Newton 修正步的齐次边界条件
    crack_nodes::Vector{Int}
    final_displacement::Float64 
end

"""
为整体式求解格式创建一个包含位移场(:u)和相场(:d)的统一自由度处理器。
"""
function create_monolithic_dofhandler_mem(grid)
    dh_a = Ferrite.DofHandler(grid)
    Ferrite.add!(dh_a, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2)
    Ferrite.add!(dh_a, :d, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}())
    Ferrite.close!(dh_a)
    return dh_a
end

"""
【新增】创建一个纯相场的自由度处理器 (为了方便计算 M_d)
"""
function create_d_dofhandler_mem(grid)
    dh_d = Ferrite.DofHandler(grid)
    Ferrite.add!(dh_d, :d, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}())
    Ferrite.close!(dh_d)
    return dh_d
end

"""
为 M-EM 算法创建边界条件。
"""
function create_mem_bcs(dh, grid, crack_nodes, final_displacement = 0.0)
    top = Ferrite.getfacetset(grid, "top")
    right = Ferrite.getfacetset(grid, "right")

    # ==============================================================
    # 1. ch_a: 随时间 t 动态缩放的边界条件 (用于 Predictor 和收敛判断)
    # ==============================================================
    ch_a = Ferrite.ConstraintHandler(dh)
    
    # 顶边全固定
    Ferrite.add!(ch_a, Ferrite.Dirichlet(:u, top, (x, t) -> zeros(2), [1, 2]))
    
    # 【核心修改】：右侧施加随 t 变化的位移。
    # 因为求解器中会调用 update!(ch_a, t_cur / t_max)，这里的 t 会在 0 到 1 之间变化。
    Ferrite.add!(ch_a, Ferrite.Dirichlet(:u, right, (x, t) -> t * final_displacement, 2))
    
    # 预制裂纹相场约束
    if !isempty(crack_nodes)
        # 注意：由于这是实际状态，预制裂纹处的相场应该是 1.0 (完全破坏)
        Ferrite.add!(ch_a, Ferrite.Dirichlet(:d, Set(crack_nodes), (x, t) -> 1.0))
    end
    Ferrite.close!(ch_a)
    Ferrite.update!(ch_a, 0.0) # 初始时刻 t=0.0

    # ==============================================================
    # 2. ch_zero: 齐次边界条件 (用于 Newton 校正步求解 Δx)
    # ==============================================================
    ch_zero = Ferrite.ConstraintHandler(dh)
    
    # 顶边全固定 (0.0)
    Ferrite.add!(ch_zero, Ferrite.Dirichlet(:u, top, (x, t) -> zeros(2), [1, 2]))
    
    # 右侧受控端：必须锁定为 0.0 (Newton 修正步不能改变边界的值)
    Ferrite.add!(ch_zero, Ferrite.Dirichlet(:u, right, (x, t) -> 0.0, 2))
    
    # 预制裂纹相场约束: 增量必须为 0.0
    if !isempty(crack_nodes)
        Ferrite.add!(ch_zero, Ferrite.Dirichlet(:d, Set(crack_nodes), (x, t) -> 0.0))
    end
    Ferrite.close!(ch_zero)
    Ferrite.update!(ch_zero, 0.0)

    return ch_a, ch_zero
end

"""
一站式初始化函数（M-EM 整体式版本）。
"""
function setup_l_tension_mem(;
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8, 
)
    grid = create_grid(msh_file)
    crack_nodes = Int[] 
    
    # 1. 创建整体式 DofHandler 和纯相场 DofHandler
    dh = create_monolithic_dofhandler_mem(grid)
    dh_d = create_d_dofhandler_mem(grid)
    
    # 2. 获取 M-EM 专用的两套 ConstraintHandler
    ch_a, ch_zero = create_mem_bcs(dh, grid, crack_nodes, final_displacement)
    
    return MonolithicTensionSetup_MEM(grid, dh, dh_d, ch_a, ch_zero, crack_nodes, final_displacement)
end