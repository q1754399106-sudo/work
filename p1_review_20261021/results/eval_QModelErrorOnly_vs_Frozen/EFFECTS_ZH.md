# QModelErrorOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 126 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 2.59%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|SoftLowQ|0.950922|0.914315|3.85|1.299817|5/12|
|LowerSNR|ThinHighFrequency|0.789360|0.751412|4.81|0.767113|0/12|
|LowerSNR|VelocityReversal|0.668010|0.629904|5.70|0.646931|0/12|
|Main|SoftLowQ|0.773074|0.760040|1.69|1.278946|6/30|
|Main|ThinHighFrequency|0.729228|0.708676|2.82|0.741977|0/30|
|Main|VelocityReversal|0.555080|0.535376|3.55|0.561290|0/30|

未通过项目：
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.6680308278224457, "relative_improvement": 0.025901230153030363}

触边/目标失配诊断：{'residual target mismatch': 65}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
