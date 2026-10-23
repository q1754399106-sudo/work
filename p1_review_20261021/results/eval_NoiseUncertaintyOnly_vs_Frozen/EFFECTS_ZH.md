# NoiseUncertaintyOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 126 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 4.75%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|SoftLowQ|0.950922|0.806014|15.24|1.087617|1/12|
|LowerSNR|ThinHighFrequency|0.789360|0.781309|1.02|0.794125|0/12|
|LowerSNR|VelocityReversal|0.668010|0.656667|1.70|0.673364|0/12|
|Main|SoftLowQ|0.773074|0.690875|10.63|0.871805|0/30|
|Main|ThinHighFrequency|0.729228|0.721803|1.02|0.748537|0/30|
|Main|VelocityReversal|0.555080|0.547021|1.45|0.566655|0/30|

未通过项目：
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.6532331557906287, "relative_improvement": 0.047478668681994574}

触边/目标失配诊断：{'residual target mismatch': 46}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
