using PhaseFieldFracture
using Ferrite
using FerriteGmsh

# 快速示例：从项目根目录运行
#   julia --project=. scripts/run_simple_tension.jl
"""
setup = setup_square_tension(
    cells = (100, 100),
    final_displacement = 0.01, # 总拉伸 0.01 mm
    crack_y = 0.0,           # 裂纹高度在正中间 y = 0.0
    crack_x_min = -1.0,       # 裂纹从左边界 x = -1.0 开始
    crack_x_max = 0.0        # 一直切到中心 x = 0.0
)
"""
setup = setup_l_tension(
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8, # 总拉伸 0.8 mm
)

# 把网格写出到 VTK 文件，方便用 ParaView 可视化检查网格质量和边界条件设置。
mkpath("data/mesh")
let dh = DofHandler(setup.grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    VTKGridFile("data/mesh/my_mesh", dh) do vtk
        # 空闭包 — 只导出网格，不写解字段
    end
end

mat = PhaseFieldMaterial(
    E = 25840, # 杨氏模量，单位 MPa
    ν = 0.18,
    gc = 0.65, # 断裂能 gc 的单位是 N/mm
    l = 10, # 尺度参数 l 建议不小于两倍网格尺寸 h
    k_tol = 1e-8
)

disp, force = solve_staggered(setup, mat; n_steps = 100, max_iter = 1000)

# 可视化载荷-位移曲线
using CairoMakie

fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1],
    xlabel = "Displacement (mm)",
    ylabel = "Reaction Force (N)",
    title = "L-shape Tension — Load‑Displacement Curve",
)

# 因为位移是向下的，也就是负值。
lines!(ax, -disp, -force; linewidth = 2, color = :steelblue)

mkpath("data/plots")
save("data/plots/load_displacement.png", fig)

println("大功告成！载荷-位移曲线已保存至 data/plots/load_displacement.png。")
println("你也可以用 ParaView 打开 data/sims/ 下的 .vtu 文件看裂纹扩展了！")
