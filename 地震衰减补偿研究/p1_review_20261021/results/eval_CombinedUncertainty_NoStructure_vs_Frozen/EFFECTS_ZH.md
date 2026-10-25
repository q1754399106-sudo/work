# CombinedUncertainty_NoStructure 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 126 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 0.01%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|SoftLowQ|0.950922|1.257294|-32.22|1.885704|8/12|
|LowerSNR|ThinHighFrequency|0.789360|0.761904|3.48|0.778585|0/12|
|LowerSNR|VelocityReversal|0.668010|0.653758|2.13|0.672657|0/12|
|Main|SoftLowQ|0.773074|0.809121|-4.66|1.493711|6/30|
|Main|ThinHighFrequency|0.729228|0.708307|2.87|0.738431|0/30|
|Main|VelocityReversal|0.555080|0.539784|2.76|0.559926|0/30|

未通过项目：
- LowerSNR/SoftLowQ/mean_noninferiority：{"Suite": "LowerSNR", "Model": "SoftLowQ", "n": 12, "baseline": 0.9509218869184454, "candidate": 1.2572937206514054, "improvement_percent": -32.2184017370541, "baseline_median": 1.0027183788149356, "candidate_median": 1.37157351033829, "baseline_p90": 1.110441033896145, "candidate_p90": 1.685850264540663, "baseline_p95": 1.1850384365551052, "candidate_p95": 1.7780010579419774, "baseline_worst": 1.27508821645001, "candidate_worst": 1.88570394131872, "baseline_failures": 6, "candidate_failures": 8}
- LowerSNR/SoftLowQ/failure_count：{"baseline": 6, "candidate": 8, "n": 12}
- Main/SoftLowQ/mean_noninferiority：{"Suite": "Main", "Model": "SoftLowQ", "n": 30, "baseline": 0.7730739096089934, "candidate": 0.8091207120891872, "improvement_percent": -4.6627886457099255, "baseline_median": 0.670037211809456, "candidate_median": 0.6645035811486879, "baseline_p90": 1.1772965998966671, "candidate_p90": 1.449841725539663, "baseline_p95": 1.2349091259273033, "candidate_p95": 1.4820323943452558, "baseline_worst": 1.2702807129234, "candidate_worst": 1.49371094670029, "baseline_failures": 6, "candidate_failures": 6}
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.6857372631678202, "relative_improvement": 8.233651824673949e-05}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
