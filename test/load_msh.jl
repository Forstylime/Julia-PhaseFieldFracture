using Ferrite
using FerriteGmsh

grid = togrid("data/mesh/l_shape.msh")

## 这个测试脚本的目的是检查从 Gmsh 导入的网格中，节点集、侧面集（facesets）、单元集和顶点集是否正确加载，并且可以通过简单的可视化来验证边界条件定义是否正确。
# 1. 检查节点集
println("1. 节点集 (Nodesets):")
for (name, set) in grid.nodesets
    println("   - 集合名称: \"$name\", 包含 $(length(set)) 个节点")
end

# 2. 检查侧面集（原 facesets）
println("\n2. 侧面集 (Facetsets - 用于施加边界条件):")
for (name, set) in grid.facetsets
    println("   - 集合名称: \"$name\", 包含 $(length(set)) 个侧面 (FacetIndex)")
end

# 3. 检查单元集
println("\n3. 单元集 (Cellsets):")
for (name, set) in grid.cellsets
    println("   - 集合名称: \"$name\", 包含 $(length(set)) 个单元")
end

# 4. 检查顶点集（新版 API 引入，通常用于最外角点约束）
println("\n4. 顶点集 (Vertexsets):")
for (name, set) in grid.vertexsets
    println("   - 集合名称: \"$name\", 包含 $(length(set)) 个顶点")
end

println("====================================================================")
## 可视化检查边界条件定义是否正确
using FerriteViz
using GLMakie

# 假设我们在 Gmsh 里定义了一个叫 "fixed_support" 的面集和一个叫 "traction" 的面集
# 我们可以分别绘制它们进行检查

# 1. 创建画布
fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data, title = "边界条件检查 (黄色为高亮边界)")

# 2. 绘制整个网格背景（灰色）
FerriteViz.plot!(ax, grid, color = :lightgrey)

# 3. 叠加高亮你想检查的 faceset（例如 "fixed_support"，显示为黄色/橘色）
if haskey(grid.facesets, "fixed_support")
    FerriteViz.plot!(ax, grid, faceset = grid.facesets["fixed_support"], color = :orange)
else
    @warn "未找到名为 'fixed_support' 的 faceset"
end

display(fig)

##
# ===== 力–位移曲线可视化 =====

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1],
    xlabel = "Displacement (mm)",
    ylabel = "Reaction Force (N)",
    title = "Force-Displacement Curve"
)

lines!(ax, disp, force, color = :blue, linewidth = 2)

# 标注峰值力
peak_idx = argmax(force)
scatter!(ax, [disp[peak_idx]], [force[peak_idx]], color = :red, markersize = 12)
text!(ax, "($(round(disp[peak_idx], digits=4)), $(round(force[peak_idx], digits=2)))",
    position = (disp[peak_idx] + 0.01, force[peak_idx]),
    fontsize = 11, color = :red)

mkpath("data/figures")
save("data/figures/force_displacement.png", fig)
save("data/figures/force_displacement.pdf", fig)

display(fig)
println("力–位移曲线已保存到 data/figures/ 文件夹。")