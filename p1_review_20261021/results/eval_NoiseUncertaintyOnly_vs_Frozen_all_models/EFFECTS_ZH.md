# NoiseUncertaintyOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 246 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 4.75%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|DevMain|Dev_DeepContrast|0.555567|0.541940|2.45|0.579887|0/15|
|DevMain|Dev_ModerateQ|0.623227|0.609744|2.16|0.629688|0/15|
|LowerSNR|SoftLowQ|0.950922|0.806014|15.24|1.087617|1/12|
|LowerSNR|ThinHighFrequency|0.789360|0.781309|1.02|0.794125|0/12|
|LowerSNR|VelocityReversal|0.668010|0.656667|1.70|0.673364|0/12|
|Main|SoftLowQ|0.773074|0.690875|10.63|0.871805|0/30|
|Main|ThinHighFrequency|0.729228|0.721803|1.02|0.748537|0/30|
|Main|VelocityReversal|0.555080|0.547021|1.45|0.566655|0/30|
|MultiCalMain|SoftLowQ|0.895024|0.749302|16.28|0.916290|0/10|
|MultiCalMain|ThinHighFrequency|0.736741|0.730707|0.82|0.762828|0/10|
|MultiCalMain|VelocityReversal|0.550469|0.541309|1.66|0.558765|0/10|
|NewGeologyMain|Holdout_DeepLowQ|0.608064|0.592770|2.52|0.752534|0/15|
|NewGeologyMain|Holdout_HighContrastReversal|0.585202|0.576503|1.49|0.589092|0/15|
|NewGeologyMain|Holdout_StrongAtten|0.884106|0.741272|16.16|1.055602|1/15|
|NewGeologyMain|Holdout_ThinInterbed|0.624575|0.613445|1.78|0.641817|0/15|

未通过项目：
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.6532331557906287, "relative_improvement": 0.047478668681994574}

触边/目标失配诊断：{'residual target mismatch': 110}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
