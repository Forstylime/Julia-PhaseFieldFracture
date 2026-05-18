using PhaseFieldFracture

# 快速示例：从项目根目录运行
#   julia --project=. scripts/run_simple_tension.jl

setup = setup_square_tension(
    cells = (40, 40),
    top_displacement = 0.02, # 总拉伸 0.02 mm
    crack_y = 0.5,           # 裂纹高度在正中间 0.5
    crack_x_min = 0.0,       # 裂纹从左边界 x=0.0 开始
    crack_x_max = 0.5        # 一直切到中心 x=0.5
)

# 把网格写出到 VTK 文件，方便用 ParaView 可视化检查网格质量和边界条件设置。
Ferrite.write_vtk("data/mesh/daa_mesh", grid)

mat = PhaseFieldMaterial(
    E = 210e3,
    ν = 0.3,
    gc = 2.7,
    l = 0.05, # 尺度参数 l 建议不小于两倍网格尺寸 h
    k_tol = 1e-8
)

disp, force = solve_staggered(setup, mat; n_steps = 100, max_u_disp = 0.02, max_iter = 100)

println("大功告成！你可以用 ParaView 打开 data/sims/ 下的 .vtu 文件看裂纹扩展了！")
