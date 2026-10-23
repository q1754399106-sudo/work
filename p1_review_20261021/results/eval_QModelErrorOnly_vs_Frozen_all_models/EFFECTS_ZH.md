# QModelErrorOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 246 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 2.59%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|DevMain|Dev_DeepContrast|0.555567|0.526693|5.20|0.569001|0/15|
|DevMain|Dev_ModerateQ|0.623227|0.600102|3.71|0.638209|0/15|
|LowerSNR|SoftLowQ|0.950922|0.914315|3.85|1.299817|5/12|
|LowerSNR|ThinHighFrequency|0.789360|0.751412|4.81|0.767113|0/12|
|LowerSNR|VelocityReversal|0.668010|0.629904|5.70|0.646931|0/12|
|Main|SoftLowQ|0.773074|0.760040|1.69|1.278946|6/30|
|Main|ThinHighFrequency|0.729228|0.708676|2.82|0.741977|0/30|
|Main|VelocityReversal|0.555080|0.535376|3.55|0.561290|0/30|
|MultiCalMain|SoftLowQ|0.895024|0.844326|5.66|1.505152|2/10|
|MultiCalMain|ThinHighFrequency|0.736741|0.716846|2.70|0.756193|0/10|
|MultiCalMain|VelocityReversal|0.550469|0.531529|3.44|0.555776|0/10|
|NewGeologyMain|Holdout_DeepLowQ|0.608064|0.600677|1.21|0.747855|0/15|
|NewGeologyMain|Holdout_HighContrastReversal|0.585202|0.561606|4.03|0.576802|0/15|
|NewGeologyMain|Holdout_StrongAtten|0.884106|0.866064|2.04|1.858446|3/15|
|NewGeologyMain|Holdout_ThinInterbed|0.624575|0.597175|4.39|0.631847|0/15|

未通过项目：
- NewGeologyMain/Holdout_ThinInterbed/WeakEventError：{"baseline": 1.4179925427020568, "candidate": 1.5124526792886868, "allowance": 0.07089962713510284, "n": 3}
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.6680308278224457, "relative_improvement": 0.025901230153030363}

触边/目标失配诊断：{'residual target mismatch': 144}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
