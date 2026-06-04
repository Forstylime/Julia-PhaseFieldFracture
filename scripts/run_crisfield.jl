using PhaseFieldFracture
using Ferrite
using FerriteGmsh
using CairoMakie

# 从项目根目录运行:
#   julia --project=. scripts/run_crisfield.jl

setup = setup_l_tension_monolithic(
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8,
)

mat = PhaseFieldMaterial(
    E = 25840,
    ν = 0.18,
    gc = 0.65,
    l = 10,
    k_tol = 1e-6,
)

disp, force, psi_energy, gf_energy = solve_crisfield(
    setup, mat;
    ρ_init = 0.02,        # 初始弧长步长 (更小的步长以解析峰值)
    max_steps = 200,      # 足够的步数以解析峰值和软化段
    tol = 1e-5,           # RMS残差容差
    max_newton = 20,      # Newton 迭代上限 (非线性区域需要更多迭代)
    λ_max = 1.0,          # 加载到最终位移
    output_freq = 5,
)

# 可视化载荷-位移曲线和能量演变
disp_plot = -disp
force_plot = force

peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷 (Crisfield): F_max = $(round(peak_force, digits=3)) N @ ū = $(round(peak_disp, digits=4)) mm")

mkpath("data/plots")

# 载荷-位移曲线
fig_load = Figure(size = (600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title = "Load-Displacement (Crisfield)",
    limits = ((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth = 2, color = :red, linestyle = :solid)
save("data/plots/load_displacement_crisfield.png", fig_load)

# 能量演变
fig_energy = Figure(size = (600, 400))
ax_energy = Axis(fig_energy[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = "Energy [N·mm]",
    title = "Energy Evolution (Crisfield)",
    limits = ((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
#lines!(ax_energy, disp_plot, psi_energy; linewidth = 2, color = :steelblue, label = L"\Psi\ \mathrm{(elastic)}")
lines!(ax_energy, disp_plot, gf_energy; linewidth = 2, color = :darkorange, label = L"\mathcal{G}_f\ \mathrm{(surface)}")
axislegend(ax_energy; position = :lt)
save("data/plots/energy_evolution_crisfield.png", fig_energy)

println("载荷-位移曲线已保存至 data/plots/load_displacement_crisfield.png。")
println("能量演变曲线已保存至 data/plots/energy_evolution_crisfield.png。")

using JLD2
@save "data/jld2/crisfield_results.jld2" disp force psi_energy gf_energy
println("Crisfield 仿真数据已保存至 data/jld2/crisfield_results.jld2。")
