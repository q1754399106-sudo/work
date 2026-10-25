# 保守策略试验结果

**结论：不提升为推荐默认版。** +3% 速度条件下可以避免严重失控，但 −3% 条件下过度拒绝了仍有用的补偿结果。

状态 COMPLETED；耗时 283.5 秒；140 行指标。

GuardedNoise 实际补偿 10/14 组观测，回退 4/14 组。

所有分支残差或标志失败 0 行；执行后核心源码哈希变化 0 个文件。

## 主场景平均误差

| 观测种子 | Frozen | 仅噪声 | 联合 | 保守策略 | 未补偿 | 保守策略相对 Frozen 改善 |
|---|---:|---:|---:|---:|---:|---:|
| 20261101 | 0.750662 | 0.673960 | 0.674592 | 0.675070 | 1.008100 | +10.07% |
| 20261107 | 0.737284 | 0.664070 | 0.662649 | 0.665271 | 1.008178 | +9.77% |

## 速度压力场景

| 条件 | 种子 | Frozen | 仅噪声 | 联合 | 保守策略 | 回退 |
|---|---|---:|---:|---:|---:|---:|
| StressVelocityMinus3 | 20261101 | 0.757030 | 0.749411 | 0.730464 | 1.000817 | 1 |
| StressVelocityMinus3 | 20261107 | 0.744276 | 0.735742 | 0.715399 | 1.000973 | 1 |
| StressVelocityPlus3 | 20261101 | 7.471158 | 7.643502 | 8.497215 | 1.000817 | 1 |
| StressVelocityPlus3 | 20261107 | 8.543862 | 8.568109 | 9.777630 | 1.000973 | 1 |

## 局部指标与校准诊断

- 种子 20261101 / CombinedUncertainty：弱反射误差 0.417778，断层邻域误差 0.673120，AVO 梯度误差 0.633707。
- 种子 20261101 / Frozen：弱反射误差 0.430321，断层邻域误差 0.669835，AVO 梯度误差 0.630645。
- 种子 20261101 / GuardedNoise：弱反射误差 0.430321，断层邻域误差 0.669835，AVO 梯度误差 0.647509。
- 种子 20261101 / NoiseUncertaintyOnly：弱反射误差 0.423477，断层邻域误差 0.670633，AVO 梯度误差 0.647509。
- 种子 20261107 / CombinedUncertainty：弱反射误差 0.461062，断层邻域误差 0.669657，AVO 梯度误差 0.613773。
- 种子 20261107 / Frozen：弱反射误差 0.471932，断层邻域误差 0.670752，AVO 梯度误差 0.614559。
- 种子 20261107 / GuardedNoise：弱反射误差 0.471932，断层邻域误差 0.670752，AVO 梯度误差 0.629171。
- 种子 20261107 / NoiseUncertaintyOnly：弱反射误差 0.467871，断层邻域误差 0.668443，AVO 梯度误差 0.629171。
- velocity_stress_StressVelocityMinus3_SoftLowQ_20261030.json：速度系数 0.97，Q=32.000000；可靠性：{"reliable": false, "reasons": ["frequency_search_boundary", "optimizer_not_converged", "calibration_fit_mismatch"], "fit_deviance": 8.024075891233855, "fit_deviance_max": 3, "interpretation": "Engineering rejection rule; selected after previous exploratory results, not a calibrated statistical test or blind evaluation"}。
- velocity_stress_StressVelocityPlus3_SoftLowQ_20261030.json：速度系数 1.03，Q=15.014928；可靠性：{"reliable": false, "reasons": ["Q_search_boundary", "optimizer_not_converged", "calibration_fit_mismatch"], "fit_deviance": 8.014623191746043, "fit_deviance_max": 3, "interpretation": "Engineering rejection rule; selected after previous exploratory results, not a calibrated statistical test or blind evaluation"}。

## 验收与补偿覆盖

- acceptance_guarded：FAIL。
  - StressVelocityMinus3/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityMinus3", "Model": "SoftLowQ", "n": 2, "baseline": 0.75065288041706, "candidate": 1.0008946401813752, "improvement_percent": -33.33654826253139, "baseline_median": 0.75065288041706, "candidate_median": 1.0008946401813752, "baseline_p90": 0.7557543590925536, "candidate_p90": 1.000957057503651, "baseline_p95": 0.7563920439269903, "candidate_p95": 1.0009648596689356, "baseline_worst": 0.757029728761427, "candidate_worst": 1.00097266183422, "baseline_failures": 0, "candidate_failures": 2}
  - StressVelocityMinus3/SoftLowQ/failure_count：{"baseline": 0, "candidate": 2, "n": 2}
  - StressVelocityMinus3/SoftLowQ/EventAmplitudeError：{"baseline": 0.709522710358759, "candidate": 0.997668500515875, "allowance": 0.035476135517937955, "n": 2}
  - StressVelocityMinus3/SoftLowQ/AVOInterceptError：{"baseline": 0.6872197787080555, "candidate": 0.991792876406848, "allowance": 0.03436098893540278, "n": 2}
  - StressVelocityMinus3/SoftLowQ/AVOGradientError：{"baseline": 0.719473474394835, "candidate": 0.997829253072292, "allowance": 0.03597367371974175, "n": 2}
  - candidate_compensation_coverage：{"compensated_observations": 10, "total": 14, "note": "Pass-through rejection is not successful attenuation compensation."}

## 解释边界

两个相反速度条件的校准归一化损失都约为 8，但直接联合反演一个严重失控、另一个仍有用。因此该损失不能单独充当重建质量判据；不能通过调一个门槛就宣称解决问题。

异常校准的回退只是原始观测透传，没有恢复衰减损失。必须把回退率与补偿误差一起报告。

Section 数据保留 Frozen，局部指标相同来自保守分支选择，不是获得新的弱反射恢复能力。

只测一个模型、一次校准、两个观测种子与两个速度条件；不提供普遍鲁棒性或统计显著性结论。

校准门槛是参考此前失败制定的探索性工程规则，尚需在独立模型与新校准上验证误拒率、漏拒率。
