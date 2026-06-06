# 算法运行结果记录

## Standard Staggered (完成)
仿真结束！VTK 文件保存在 data/sims/staggered 目录下。
总Newton迭代次数: 4708
计算耗时: 285.17 秒
峰值载荷: F_max = 328.4289 N @ ū = 0.696 mm
载荷-位移曲线已保存至 data/plots/load_displacement_staggered.png。
能量演变曲线已保存至 data/plots/energy_evolution_staggered.png。
Staggered 仿真数据已保存至 data/jld2/staggered_results.jld2。

## S-EM (完成)
S-EM 仿真结束！VTK 文件保存在 data/sims/sem 目录下。
总载荷步数: 169
总Newton迭代次数: 9787
计算耗时: 1922.38 秒
峰值载荷 (S-EM): F_max = 328.341 N @ ū = 0.6935 mm
载荷-位移曲线已保存至 data/plots/load_displacement_sem.png。
能量演变曲线已保存至 data/plots/energy_evolution_sem.png。
S-EM 仿真数据已保存至 data/jld2/sem_results.jld2。

## Crisfield Arc Length (已修复 ✓)
=== 载荷步 265 | λ = 0.9979 | ρ = 0.05 ===
仿真结束！VTK 文件保存在 data/sims/crisfield 目录下。
计算耗时: 82.11 秒
峰值载荷 (Crisfield): F_max = 328.412 N @ ū = 0.6971 mm
载荷-位移曲线已保存至 data/plots/load_displacement_crisfield.png。
能量演变曲线已保存至 data/plots/energy_evolution_crisfield.png。
Crisfield 仿真数据已保存至 data/jld2/crisfield_results.jld2。
