using PhaseFieldFracture
using Ferrite
using FerriteGmsh
using CairoMakie

# 从项目根目录运行:
#   julia --project=. scripts/run_crisfield.jl

setup = setup_tension_monolithic(
    msh_file = "data/mesh/square.msh",
    final_displacement = 0.3,
    dir = 1,
    fixed_face = "left"
)
"""
mat = PhaseFieldMaterial(
    E = 25840,
    ν = 0.18,
    gc = 0.65,
    l = 10,
    k_tol = 1e-8
)
"""
mat = PhaseFieldMaterial(
    E = 100,
    ν = 0.3,
    gc = 1.0,
    l = 0.05,
    k_tol = 1e-8
)

disp, force, psi_energy, gf_energy = solve_crisfield(
    setup, mat;
    ρ_init = 0.01,
    max_steps = 400,
    tol = 1e-5,
    max_newton = 20,
    output_freq = 10,
    λ_max = 1.0,
    enforce_irreversibility = true
)

# 可视化载荷-位移曲线和能量演变
disp_plot = disp
force_plot = force

peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷 (Crisfield): F_max = $(round(peak_force, digits=3)) N @ ū = $(round(peak_disp, digits=4)) mm")

mkpath("data/plots2")

# 载荷-位移曲线
fig_load = Figure(size = (600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title = "Load-Displacement (Crisfield)",
    #limits = ((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth = 2, color = :red, linestyle = :solid)
save("data/plots2/load_displacement_crisfield.png", fig_load)

# 能量演变
fig_energy = Figure(size = (600, 400))
ax_energy = Axis(fig_energy[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = "Energy [N·mm]",
    title = "Energy Evolution (Crisfield)",
    #limits = ((0, maximum(disp_plot)), (0, 82)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
#lines!(ax_energy, disp_plot, psi_energy; linewidth = 2, color = :steelblue, label = L"\Psi\ \mathrm{(elastic)}")
lines!(ax_energy, disp_plot, gf_energy; linewidth = 2, color = :darkorange, label = L"\mathcal{G}_f\ \mathrm{(surface)}")
axislegend(ax_energy; position = :lt)
save("data/plots2/energy_evolution_crisfield.png", fig_energy)

println("载荷-位移曲线已保存至 data/plots2/load_displacement_crisfield.png。")
println("能量演变曲线已保存至 data/plots2/energy_evolution_crisfield.png。")

using JLD2
@save "data/jld2/crisfield_results2.jld2" disp force psi_energy gf_energy
println("Crisfield 仿真数据已保存至 data/jld2/crisfield_results2.jld2。")
