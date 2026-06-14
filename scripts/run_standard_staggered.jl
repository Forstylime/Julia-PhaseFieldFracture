using PhaseFieldFracture
using Ferrite
using FerriteGmsh

# 快速示例：从项目根目录运行
#   julia --project=. scripts/run_standard_staggered.jl

"""
setup = setup_l_tension(
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = -0.8, # 总拉伸 -0.8 mm (向下)
)
    
mat = PhaseFieldMaterial(
    E = 25840, # 杨氏模量，单位 MPa
    ν = 0.18,
    gc = 0.65, # 断裂能 gc 的单位是 N/mm
    l = 10, # 尺度参数 l 建议不小于两倍网格尺寸 h
    k_tol = 1e-8
)
"""
setup = setup_square_tension(
    msh_file = "data/mesh/square.msh",
    final_displacement = 0.3, # 总拉伸 0.3 mm (向右)
)

mat = PhaseFieldMaterial(
    E = 100,  # 杨氏模量，单位 MPa
    ν = 0.3,
    gc = 1.0,   # 断裂能 gc 的单位是 N/mm
    l = 0.05, # 尺度参数 l 建议不小于两倍网格尺寸 h
    k_tol = 1e-8
)
# 把网格写出到 VTK 文件，方便用 ParaView 可视化检查网格质量和边界条件设置。
mkpath("data/mesh")
let dh = DofHandler(setup.grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    VTKGridFile("data/mesh/square_mesh", dh) do vtk
        # 空闭包 — 只导出网格，不写解字段
    end
end

disp, force, psi_energy, gf_energy = solve_staggered(setup, mat; 
    n_steps = 100, 
    max_iter = 2000, 
    enforce_irreversibility = false
)

# 可视化载荷-位移曲线 and 能量演变
using CairoMakie

disp_plot = disp
force_plot = force

# 找到峰值载荷及其对应位移
peak_idx = argmax(force_plot)
peak_force = force_plot[peak_idx]
peak_disp = disp_plot[peak_idx]
println("峰值载荷: F_max = $(round(peak_force, digits=4)) N @ ū = $(round(peak_disp, digits=4)) mm")

# 图一：载荷-位移曲线
fig_load = Figure(size = (600, 400))
ax_load = Axis(fig_load[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = L"F_{\mathrm{reaction}}~\mathrm{[N]}",
    title = "Load‑Displacement",
    #limits = ((0, maximum(disp_plot)), (0, nothing)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
lines!(ax_load, disp_plot, force_plot; linewidth = 2, color = :red, linestyle = :dash) # 红色虚线
save("data/plots2/load_displacement_staggered.png", fig_load)

# 图二：能量演变
fig_energy = Figure(size = (600, 400))
ax_energy = Axis(fig_energy[1, 1],
    xlabel = L"\bar{u}~\mathrm{[mm]}",
    ylabel = "Energy [N·mm]",
    title = "Energy Evolution",
    limits = ((0, maximum(disp_plot)), (0, 1.5)),
    xgridvisible = true,
    ygridvisible = true,
    xgridcolor = :lightgray,
    ygridcolor = :lightgray,
)
#lines!(ax_energy, disp_plot, psi_energy; linewidth = 2, color = :steelblue, label = L"\Psi\ \mathrm{(elastic)}")
lines!(ax_energy, disp_plot, gf_energy; linewidth = 2, color = :darkorange, label = L"\mathcal{G}_f\ \mathrm{(surface)}")
axislegend(ax_energy; position = :lt)
save("data/plots2/energy_evolution_staggered.png", fig_energy)

println("载荷-位移曲线已保存至 data/plots/load_displacement_staggered.png。")
println("能量演变曲线已保存至 data/plots/energy_evolution_staggered.png。")

# 保存数据以供后续分析
using JLD2
@save "data/jld2/staggered_results2.jld2" disp force psi_energy gf_energy
println("Staggered 仿真数据已保存至 data/jld2/staggered_results2.jld2。")