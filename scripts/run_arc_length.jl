using PhaseFieldFracture
using Ferrite
using FerriteGmsh
using CairoMakie

println("Setting up L-shape test...")
# L-shape test case, standard values
mat = PhaseFieldFracture.PhaseFieldMaterial(
    E=25840, # MPa
    ν=0.18,
    gc=0.65, # kN/mm
    l=10,   # mm
    k_tol=1e-6
)

# We will use the setup_l_tension to get grid, boundary sets, etc.
# The setup function gives final_displacement = 0.3 or 1.0 depending on t_max
# We'll set final_displacement = 1.0 and run until t_max = 1.0.
setup = setup_l_tension(
    msh_file="data/mesh/l_shape.msh",
    final_displacement=0.8
)

# Run the arc-length solver
lambdas, displacements, reaction_forces, psi_energy, gf_energy = solve_arc_length(
    setup, mat,
    n_steps_max=100,
    ρ_init=0.04,
    tol_newton=1e-3,
    max_newton_iter=50,
    t_max=1.0,
    output_freq=10
)

println("Load-Displacement Curve Data:")
for i in 1:length(lambdas)
    println("λ = $(round(lambdas[i], digits=4)), u = $(round(displacements[i], digits=4)), f = $(round(reaction_forces[i], digits=4))")
end
# 可视化载荷-位移曲线和能量演变
disp_plot = displacements
force_plot = reaction_forces

peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷 (S-EM): F_max = $(round(peak_force, digits=3)) N @ ū = $(round(peak_disp, digits=4)) mm")

mkpath("data/plots")

# 图一：载荷-位移曲线
fig_load = Figure(size=(600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel=L"\bar{u}~\mathrm{[mm]}",
    ylabel=L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title="Load-Displacement (S-EM)",
    limits=((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible=true,
    ygridvisible=true,
    xgridcolor=:lightgray,
    ygridcolor=:lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth=2, color=:red, linestyle=:solid) # 红色实线
save("data/plots/load_displacement_crisfield.png", fig_load)

# 图二：能量演变
fig_energy = Figure(size=(600, 400))
ax_energy = Axis(fig_energy[1, 1],
    xlabel=L"\bar{u}~\mathrm{[mm]}",
    ylabel="Energy [N·mm]",
    title="Energy Evolution (S-EM)",
    limits=((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible=true,
    ygridvisible=true,
    xgridcolor=:lightgray,
    ygridcolor=:lightgray,
)
lines!(ax_energy, disp_plot, psi_energy; linewidth=2, color=:steelblue, label=L"\Psi\ \mathrm{(elastic)}")
lines!(ax_energy, disp_plot, gf_energy; linewidth=2, color=:darkorange, label=L"\mathcal{G}_f\ \mathrm{(surface)}")
axislegend(ax_energy; position=:lt)
save("data/plots/energy_evolution_crisfield.png", fig_energy)

println("载荷-位移曲线已保存至 data/plots/load_displacement_crisfield.png。")
println("能量演变曲线已保存至 data/plots/energy_evolution_crisfield.png。")

println("Done!")

