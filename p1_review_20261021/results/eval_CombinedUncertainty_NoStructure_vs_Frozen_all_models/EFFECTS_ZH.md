# CombinedUncertainty_NoStructure 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 246 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 0.01%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|DevMain|Dev_DeepContrast|0.555567|0.531427|4.35|0.574304|0/15|
|DevMain|Dev_ModerateQ|0.623227|0.599749|3.77|0.617285|0/15|
|LowerSNR|SoftLowQ|0.950922|1.257294|-32.22|1.885704|8/12|
|LowerSNR|ThinHighFrequency|0.789360|0.761904|3.48|0.778585|0/12|
|LowerSNR|VelocityReversal|0.668010|0.653758|2.13|0.672657|0/12|
|Main|SoftLowQ|0.773074|0.809121|-4.66|1.493711|6/30|
|Main|ThinHighFrequency|0.729228|0.708307|2.87|0.738431|0/30|
|Main|VelocityReversal|0.555080|0.539784|2.76|0.559926|0/30|
|MultiCalMain|SoftLowQ|0.895024|0.947290|-5.84|1.763641|2/10|
|MultiCalMain|ThinHighFrequency|0.736741|0.716480|2.75|0.753852|0/10|
|MultiCalMain|VelocityReversal|0.550469|0.536720|2.50|0.554141|0/10|
|NewGeologyMain|Holdout_DeepLowQ|0.608064|0.591170|2.78|0.757992|0/15|
|NewGeologyMain|Holdout_HighContrastReversal|0.585202|0.567343|3.05|0.579093|0/15|
|NewGeologyMain|Holdout_StrongAtten|0.884106|0.934806|-5.73|1.970200|3/15|
|NewGeologyMain|Holdout_ThinInterbed|0.624575|0.596875|4.44|0.623597|0/15|

未通过项目：
- LowerSNR/SoftLowQ/mean_noninferiority：{"Suite": "LowerSNR", "Model": "SoftLowQ", "n": 12, "baseline": 0.9509218869184454, "candidate": 1.2572937206514054, "improvement_percent": -32.2184017370541, "baseline_median": 1.0027183788149356, "candidate_median": 1.37157351033829, "baseline_p90": 1.110441033896145, "candidate_p90": 1.685850264540663, "baseline_p95": 1.1850384365551052, "candidate_p95": 1.7780010579419774, "baseline_worst": 1.27508821645001, "candidate_worst": 1.88570394131872, "baseline_failures": 6, "candidate_failures": 8}
- LowerSNR/SoftLowQ/failure_count：{"baseline": 6, "candidate": 8, "n": 12}
- Main/SoftLowQ/mean_noninferiority：{"Suite": "Main", "Model": "SoftLowQ", "n": 30, "baseline": 0.7730739096089934, "candidate": 0.8091207120891872, "improvement_percent": -4.6627886457099255, "baseline_median": 0.670037211809456, "candidate_median": 0.6645035811486879, "baseline_p90": 1.1772965998966671, "candidate_p90": 1.449841725539663, "baseline_p95": 1.2349091259273033, "candidate_p95": 1.4820323943452558, "baseline_worst": 1.2702807129234, "candidate_worst": 1.49371094670029, "baseline_failures": 6, "candidate_failures": 6}
- MultiCalMain/SoftLowQ/mean_noninferiority：{"Suite": "MultiCalMain", "Model": "SoftLowQ", "n": 10, "baseline": 0.8950243097772571, "candidate": 0.9472900645457073, "improvement_percent": -5.839590522569993, "baseline_median": 0.7565734266471069, "candidate_median": 0.7609860977477495, "baseline_p90": 1.393400448576674, "candidate_p90": 1.63532459196835, "baseline_p95": 1.4273241193316868, "candidate_p95": 1.6994827425092647, "baseline_worst": 1.4612477900867, "candidate_worst": 1.76364089305018, "baseline_failures": 2, "candidate_failures": 2}
- MultiCalMain/VelocityReversal/WeakEventError：{"baseline": 0.5175931509621565, "candidate": 0.5551928926161485, "allowance": 0.025879657548107827, "n": 2}
- NewGeologyMain/Holdout_StrongAtten/mean_noninferiority：{"Suite": "NewGeologyMain", "Model": "Holdout_StrongAtten", "n": 15, "baseline": 0.8841062527400156, "candidate": 0.9348060530090573, "improvement_percent": -5.734582253197873, "baseline_median": 0.710040099814264, "candidate_median": 0.710096760060483, "baseline_p90": 1.564419795762366, "candidate_p90": 1.853136306570776, "baseline_p95": 1.6867045069248507, "candidate_p95": 1.912535948081557, "baseline_worst": 1.83170616586718, "candidate_worst": 1.97019998046471, "baseline_failures": 3, "candidate_failures": 3}
- main_equal_model_improvement：{"baseline": 0.6857937290357045, "candidate": 0.6857372631678202, "relative_improvement": 8.233651824673949e-05}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
