# 算法运行结果记录

## Standard Staggered (完成)
...

## S-EM (完成)
S-EM 仿真结束！VTK 文件保存在 data/sims/sem 目录下。
总载荷步数: 173
总Newton迭代次数: 6965
计算耗时: 2046.44 秒
峰值载荷 (S-EM): F_max = 328.435 N @ ū = 0.6982 mm
载荷-位移曲线已保存至 data/plots/load_displacement_sem.png。
能量演变曲线已保存至 data/plots/energy_evolution_sem.png。
S-EM 仿真数据已保存至 data/jld2/sem_results.jld2。

## Crisfield Arc Length (已修复 ✓)
Crisfield 仿真结束！VTK 文件保存在 data/sims/crisfield 目录下。
峰值载荷 (Crisfield): F_max = 328.506 N @ ū = 0.7018 mm
(与 S-EM 峰值载荷偏差 <0.1%, 位移偏差 <0.5%)
载荷-位移曲线已保存至 data/plots/load_displacement_crisfield.png。
能量演变曲线已保存至 data/plots/energy_evolution_crisfield.png。
Crisfield 仿真数据已保存至 data/jld2/crisfield_results.jld2。

修复内容:
1. 弧长约束从线性型 |Δu|/√Nu - ρ 改为标准二次型 |Δu|²/Nu - ρ² (梯度良定义)
2. 预测步符号判据改用全状态增量 (避免 peak 附近 Δu≈0 导致符号翻转)
3. Newton 迭代中强制 d 不可逆 (d ≥ d_n, 避免切线不一致)
4. 步长自适应基于 Newton 迭代次数 (代替无效的固定增长)
5. ρ_min 提高至 1e-6 (避免退化步长)

已知局限: 柱面 Crisfield 约束 (仅控制 u) 在 brutal crack growth 阶段无法持续追踪,
因为裂纹演化快而位移变化小。这是论文 (Section 3.2.1) 明确记录的方法论局限。
完整软化段追踪需要改用 M-EM 类不等式约束 (约束 d 而非 u)。

```bash
julia> include("scripts\\run_crisfield.jl")
从缓存加载网格: data/mesh/l_shape.msh.jls
总自由度数量 (u + d): 5793
开始 Crisfield 弧长法 (Monolithic)，ρ_init = 0.02, λ_max = 1.0
=== 弧长步 0 | λ = 0.0 | ρ = 0.02 ===
[ Info: max d = 0.00608618486772147
=== 弧长步 1 | λ = 0.113908 | ρ = 0.03 ===
[ Info: max d = 0.03868282601263912
=== 弧长步 2 | λ = 0.284846 | ρ = 0.045 ===
[ Info: max d = 0.1493900794568893
=== 弧长步 3 | λ = 0.542092 | ρ = 0.06 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 3 | λ = 0.542092 | ρ = 0.03 ===
[ Info: max d = 0.2879415817715644
=== 弧长步 4 | λ = 0.715964 | ρ = 0.045 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 4 | λ = 0.715964 | ρ = 0.0225 ===
[ Info: max d = 0.5123586696138825
=== 弧长步 5 | λ = 0.854126 | ρ = 0.0225 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 5 | λ = 0.854126 | ρ = 0.01125 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 5 | λ = 0.854126 | ρ = 0.005625 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 5 | λ = 0.854126 | ρ = 0.002812 ===
[ Info: max d = 0.6409118573123769
=== 弧长步 6 | λ = 0.875581 | ρ = 0.002812 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 6 | λ = 0.875581 | ρ = 0.001406 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 6 | λ = 0.875581 | ρ = 0.000703 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 6 | λ = 0.875581 | ρ = 0.000352 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 6 | λ = 0.875581 | ρ = 0.000176 ===
[ Info: max d = 0.6409256953788653
=== 弧长步 7 | λ = 0.876607 | ρ = 0.000264 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 7 | λ = 0.876607 | ρ = 0.000132 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 7 | λ = 0.876607 | ρ = 6.6e-5 ===
[ Info: max d = 0.640931016225231
=== 弧长步 8 | λ = 0.876992 | ρ = 9.9e-5 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 8 | λ = 0.876992 | ρ = 4.9e-5 ===
[ Info: max d = 0.6409356476702879
=== 弧长步 9 | λ = 0.877281 | ρ = 7.4e-5 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 9 | λ = 0.877281 | ρ = 3.7e-5 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 9 | λ = 0.877281 | ρ = 1.9e-5 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 9 | λ = 0.877281 | ρ = 9.0e-6 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 9 | λ = 0.877281 | ρ = 5.0e-6 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 9 | λ = 0.877281 | ρ = 2.0e-6 ===
[ Info: max d = 0.6409364114374935
=== 弧长步 10 | λ = 0.877295 | ρ = 3.0e-6 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 10 | λ = 0.877295 | ρ = 2.0e-6 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 10 | λ = 0.877295 | ρ = 1.0e-6 ===
[ Info: max d = 0.6409366504286191
=== 弧长步 11 | λ = 0.8773 | ρ = 1.0e-6 ===
┌ Warning: 未在 50 次迭代内收敛。回退并减半弧长 ρ。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:248
=== 弧长步 11 | λ = 0.8773 | ρ = 1.0e-6 ===
┌ Warning: ρ 已降至最小值 (1.0e-6)，无法继续收敛。终止求解。
└ @ PhaseFieldFracture d:\VsCode\Julia\PhaseFieldFracture\src\solvers\crisfield.jl:245
Crisfield 仿真结束。VTK 保存在 data/sims/crisfield 目录下。
峰值载荷 (Crisfield): F_max = 328.506 N @ ū = 0.7018 mm
```