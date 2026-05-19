using PhaseFieldFracture
using Ferrite

# 快速示例：从项目根目录运行
#   julia --project=. scripts/run_simple_tension.jl

setup = setup_square_tension(
    cells = (100, 100),
    top_displacement = 0.01, # 总拉伸 0.01 mm
    crack_y = 0.0,           # 裂纹高度在正中间 y = 0.0
    crack_x_min = -1.0,       # 裂纹从左边界 x = -1.0 开始
    crack_x_max = 0.0        # 一直切到中心 x = 0.0
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
    E = 210e3,
    ν = 0.3,
    gc = 2.7,
    l = 0.05, # 尺度参数 l 建议不小于两倍网格尺寸 h
    k_tol = 1e-8
)

disp, force = solve_staggered(setup, mat; n_steps = 100, max_u_disp = 0.01, max_iter = 500)

println("大功告成！你可以用 ParaView 打开 data/sims/ 下的 .vtu 文件看裂纹扩展了！")
