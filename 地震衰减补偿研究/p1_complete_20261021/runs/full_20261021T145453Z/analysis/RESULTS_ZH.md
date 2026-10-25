# P1 全量 MATLAB 重跑与补充验收

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run + validate
- Verification Status: VERIFIED for executed coverage and solver checks; scientific acceptance below
- 数值实现：用户此次附件中的 MATLAB 七分支；超参数未调整。原有九模型属于已见模型回归验证。

实际覆盖 **40,824 行指标、2,916 组配对观测**；覆盖完整：True。
9 模型 × 3 校准 × 6 观测种子；18 个场景/条件组合：Main 5、LowerSNR 5、四种压力各 2。
每组观测包含七方法和 Direct/Split 两路线；指标行数不是独立样本数。

**CombinedUncertainty 工程验收：FAIL，314/328；加入基准求解质量与完整覆盖要求后：FAIL。**

## Main 等模型加权改善与交叉重采样区间

| 模型组 | 模型数 | 改善 % | 95% 区间 % |
|---|---:|---:|---|
| historical | 3 | 5.436 | [2.927, 8.671] |
| development | 2 | 4.908 | [3.964, 5.868] |
| heldout_new_geology | 4 | 6.148 | [3.577, 9.539] |
| all | 9 | 5.647 | [3.842, 7.663] |

模型、校准、观测种子同时重采样；校准与观测按交叉设计抽样，跨模型/场景共享种子抽样权重，场景作为固定基准保留。三次校准、少量模型仍限制区间稳定性。

## 未通过项（全部保留）

- LowerSNR/Holdout_ThinInterbed/WeakEventError：{"baseline": 1.0292816908792286, "candidate": 1.1859410643091761, "allowance": 0.05146408454396143, "n": 18}
- Main/ThinHighFrequency/WeakEventError：{"baseline": 1.1419495008520673, "candidate": 1.2029134431845157, "allowance": 0.05709747504260337, "n": 18}
- StressVelocityError/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Dev_ModerateQ", "n": 36, "baseline": 12.236478804956636, "candidate": 12.803399145693646, "improvement_percent": -4.633034958613802, "baseline_median": 13.1079701613949, "candidate_median": 13.457355750030949, "baseline_p90": 15.78237672772405, "candidate_p90": 16.6376961284395, "baseline_p95": 16.010052777633923, "candidate_p95": 17.4024241709052, "baseline_worst": 17.3511560774854, "candidate_worst": 18.9056018216322, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Dev_ModerateQ/EventAmplitudeError：{"baseline": 3.237604398825965, "candidate": 3.4651123957071315, "allowance": 0.16188021994129825, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOInterceptError：{"baseline": 2.7282603219423742, "candidate": 2.982515134787779, "allowance": 0.13641301609711873, "n": 36}
- StressVelocityError/Dev_ModerateQ/AVOGradientError：{"baseline": 2.3894831285401654, "candidate": 2.8017961793107125, "allowance": 0.11947415642700827, "n": 36}
- StressVelocityError/Holdout_StrongAtten/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "Holdout_StrongAtten", "n": 36, "baseline": 14.457904392050219, "candidate": 15.913835640285011, "improvement_percent": -10.07014024131565, "baseline_median": 12.676617824936663, "candidate_median": 13.48624467105379, "baseline_p90": 23.9993398175421, "candidate_p90": 27.48626766753935, "baseline_p95": 26.6926891176284, "candidate_p95": 32.70219546471733, "baseline_worst": 26.6926892998688, "candidate_worst": 32.7091049732802, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/Holdout_StrongAtten/EventAmplitudeError：{"baseline": 3.7553818613177787, "candidate": 4.243865303375985, "allowance": 0.18776909306588896, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOInterceptError：{"baseline": 3.037779074454735, "candidate": 3.4748344437541947, "allowance": 0.15188895372273675, "n": 36}
- StressVelocityError/Holdout_StrongAtten/AVOGradientError：{"baseline": 2.281710106016304, "candidate": 2.817732306279983, "allowance": 0.11408550530081521, "n": 36}
- StressVelocityError/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "SoftLowQ", "n": 36, "baseline": 14.120586199193193, "candidate": 15.643920657594458, "improvement_percent": -10.788039794610672, "baseline_median": 13.015832476676675, "candidate_median": 14.163460797908051, "baseline_p90": 21.57603703907565, "candidate_p90": 24.47655443576945, "baseline_p95": 25.36939763507525, "candidate_p95": 30.8736716429768, "baseline_worst": 25.4147687872566, "candidate_worst": 30.9421007362676, "baseline_failures": 36, "candidate_failures": 36}
- StressVelocityError/SoftLowQ/EventAmplitudeError：{"baseline": 4.558573031099418, "candidate": 5.2378895151726095, "allowance": 0.2279286515549709, "n": 36}
- StressVelocityError/SoftLowQ/AVOInterceptError：{"baseline": 3.354832433157156, "candidate": 3.9260610167413685, "allowance": 0.1677416216578578, "n": 36}
- StressVelocityError/SoftLowQ/AVOGradientError：{"baseline": 2.1350922144952227, "candidate": 2.7646271897218186, "allowance": 0.10675461072476114, "n": 36}

## 所有方法求解质量

| 方法 | 行数 | 最终残差/标志失败 | 任一选参尺度失败 | 最大最终残差 |
|---|---:|---:|---:|---:|
| Baseline | 5832 | 0 | 0 | 2.965e-13 |
| Frozen | 5832 | 0 | 0 | 1.000e-07 |
| StrongerDamping | 5832 | 0 | 0 | 9.999e-08 |
| NoiseUncertaintyOnly | 5832 | 0 | 0 | 9.998e-08 |
| QModelErrorOnly | 5832 | 0 | 0 | 1.000e-07 |
| CombinedUncertainty | 5832 | 0 | 0 | 9.999e-08 |
| CombinedUncertainty_NoStructure | 5832 | 0 | 0 | 4.008e-13 |

## 执行与解释边界

- 四种压力分别为 30° 常相位旋转、公开速度 +3% 且层厚固定、远偏移最大 4ms 残余 NMO、AR(1) rho=0.55 有色噪声。速度条件包含相应走时失配，重新校准并重建算子。
- 本次速度定义遵循配置中 +3% 的设置；不同于旧报告 vp +4%、h +3.5% 的联合扰动，不能直接混为同一测试。
- 残余 NMO 驱动在 IFFT 后先截取观测长度，避免原驱动的维度赋值错误。
- 完整的采样计划和源码哈希在实验前保存；运行前后源码变化见 EXECUTION_SUMMARY.json。
- 匹配先验审计直接比较实际输出记录中的 Q、频率、尺度、噪声、lambda、H 范数及目标系数；同一次函数调用产生两路线。
- 生产入口投毒测试使用相同 D/calD，分别从正常/私有字段污染配置生成 public 后调用实际 p1_compensate，并核验非法真值字段会被拒绝。
- 校准触边或优化器非正退出记录：21/54；见 calibration_diagnostics.json。保留这些失败，未按真值筛除。
- 原有模型和已观察结果不会因重新冻结被变成独立盲测。新地质泛化仍需未参与开发的独立模型或真实数据。
- 相同整数随机种子在 MATLAB 与 Python 中不保证相同随机噪声；本次是完整 MATLAB 执行，不能凭旧 Python CSV 宣称逐值等价。
- 11 类统计谬误均检查：整体/分组反转、群体到个体外推、选择偏差、碰撞条件化、基率、均值回归、幸存者、多重搜索、分析路径、因果误用、反向因果。重点限制是已见模型选择、局部均值不能保证每个事件、多分支探索和有限模型数量；失败运行及全部样本已保留，未作结果驱动筛除。

完整逐模型结果见 group_summary.csv；所有未通过项见 acceptance_summary.json。
