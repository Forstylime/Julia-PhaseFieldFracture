using PhaseFieldFracture
using Ferrite
using FerriteGmsh
using CairoMakie

# 快速示例：从项目根目录运行
#   julia --project=. scripts/run_mem.jl

setup = setup_l_tension_mem(
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8,
)

mat = PhaseFieldMaterial(
    E = 25840,
    ν = 0.18,
    gc = 0.65,
    l = 10,
    k_tol = 1e-8,
)

disp, force, psi_energy, gf_energy = solve_mem(
    setup, mat;
    max_steps = 500,
    ρ_init = 0.8658,             # 相场演化上限 (弧长)
    tol_newton = 1e-4,      # Newton 迭代收敛容差
    tol_kkt = 1e-6,         # 状态切换判断容差
    max_newton = 25,        # 最大 Newton 迭代次数
    t_max = 86.58,          # 最大伪时间 (对应最终位移)
    output_freq = 10,
)
    
# 可视化载荷-位移曲线和能量演变
disp_plot = -disp
force_plot = -force

peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷 (M-EM): F_max = $(round(peak_force, digits=3)) N @ ū = $(round(peak_disp, digits=4)) mm")

mkpath("data/plots")

# 图一：载荷-位移曲线
fig_load = Figure(size = (600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title = "Load-Displacement (M-EM)",
    limits = ((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth = 2, color = :red, linestyle = :solid) # 红色实线
save("data/plots/load_displacement_mem.png", fig_load)

# 图二：能量演变
fig_energy = Figure(size = (600, 400))
ax_energy = Axis(fig_energy[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = "Energy [N·mm]",
    title = "Energy Evolution (M-EM)",
    limits = ((0, maximum(disp_plot)), (0, 82)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
#lines!(ax_energy, disp_plot, psi_energy; linewidth = 2, color = :steelblue, label = L"\Psi\ \mathrm{(elastic)}")
lines!(ax_energy, disp_plot, gf_energy; linewidth = 2, color = :blue, label = L"\mathcal{G}_f\ \mathrm{(surface)}")
axislegend(ax_energy; position = :lt)
save("data/plots/energy_evolution_mem.png", fig_energy)

println("载荷-位移曲线已保存至 data/plots/load_displacement_mem.png。")
println("能量演变曲线已保存至 data/plots/energy_evolution_mem.png。")

# 保存数据以供后续分析
using JLD2
@save "data/jld2/mem_results.jld2" disp force psi_energy gf_energy
println("M-EM 仿真数据已保存至 data/jld2/mem_results.jld2。")
