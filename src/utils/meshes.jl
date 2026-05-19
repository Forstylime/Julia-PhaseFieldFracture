"""
    refine_range(L, n; center=0.0, ratio=3.0)

生成 `n` 个节点在区间 `[-L/2, L/2]` 上的一维坐标，在 `center` 附近加密。

使用 sinh 函数实现节点从粗到细的平滑过渡：
- `ratio = 1.0`  → 均匀分布。
- `ratio > 1.0`  → `center` 附近加密。
- `ratio` 越大加密越强（建议 3.0 ~ 5.0）。
"""
function refine_range(L, n::Int; center::Float64 = 0.0, ratio::Float64 = 3.0)
    half = L / 2
    xs = range(-1.0, 1.0, length = n)
    stretched = sinh.(ratio .* xs) ./ sinh(ratio)
    return center .+ half .* stretched
end

"""
    refine_grid!(grid, x_coords, y_coords)

将已有 Ferrite 网格的节点坐标替换为指定的非均匀坐标。

假设 `grid` 是由 `nx × ny` 个节点组成的结构化四边形网格，
节点按列优先排列（与 `Ferrite.generate_grid(Quadrilateral, ...)` 一致）。
"""
function refine_grid!(grid, x_coords::Vector{Float64}, y_coords::Vector{Float64})
    nx = length(x_coords)
    ny = length(y_coords)

    for j in 1:ny
        for i in 1:nx
            idx = (j - 1) * nx + i
            grid.nodes[idx] = Ferrite.Node((x_coords[i], y_coords[j]))
        end
    end

    return grid
end
