using PhaseFieldFracture
using Ferrite
using FerriteGmsh
using CairoMakie

# 从项目根目录运行:
#   julia --project=. scripts/run_L2.jl

setup = setup_l_tension_monolithic(
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8,
)

mat = PhaseFieldMaterial(
    E = 25840,
    ν = 0.18,
    gc = 0.65,
    l = 10,
    k_tol = 1e-8
)

disp, force, psi_energy, gf_energy = solve_l2(
    setup, mat;
    ρ_init = 0.02,        # 初始弧长步长。建议 0.01 ~ 0.05 (更小的步长以解析峰值)
    ρ_max = 0.2,          # 断裂期间的最大弧长限制。越小，裂纹扩展阶段的步数越多、曲线越平滑。建议 0.1 ~ 0.5
    Δλ_max = 0.04,        # 弹性/未开裂期间的最大位移步长。越小，弹性直线上点越多。建议 0.02 ~ 0.05
    target_iters = 4.0,   # 自适应激进程度。设为 3.0 会更趋向于保守切分，设为 5.0 会更激进一些。建议 3.0 ~ 5.0
    max_steps = 400,      # 适当增加最大步数以确保高精度模式下能完成
)

# 可视化载荷-位移曲线和能量演变
disp_plot = -disp
force_plot = -force

peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷 (L2): F_max = $(round(peak_force, digits=3)) N @ ū = $(round(peak_disp, digits=4)) mm")

mkpath("data/plots")

# 载荷-位移曲线
fig_load = Figure(size = (600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title = "Load-Displacement (L2)",
    limits = ((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth = 2, color = :red, linestyle = :solid)
save("data/plots/load_displacement_l2.png", fig_load)

# 能量演变
fig_energy = Figure(size = (600, 400))
ax_energy = Axis(fig_energy[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = "Energy [N·mm]",
    title = "Energy Evolution (L2)",
    limits = ((0, maximum(disp_plot)), (0, 82)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
#lines!(ax_energy, disp_plot, psi_energy; linewidth = 2, color = :steelblue, label = L"\Psi\ \mathrm{(elastic)}")
lines!(ax_energy, disp_plot, gf_energy; linewidth = 2, color = :darkorange, label = L"\mathcal{G}_f\ \mathrm{(surface)}")
axislegend(ax_energy; position = :lt)
save("data/plots/energy_evolution_l2.png", fig_energy)

println("载荷-位移曲线已保存至 data/plots/load_displacement_l2.png。")
println("能量演变曲线已保存至 data/plots/energy_evolution_l2.png。")

using JLD2
@save "data/jld2/l2_results.jld2" disp force psi_energy gf_energy
println("L2 仿真数据已保存至 data/jld2/l2_results.jld2。")
