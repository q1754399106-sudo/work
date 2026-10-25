# ConfidenceBalanced 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 48 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 1.83%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.668043|2.19|0.719855|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.713048|1.60|0.728633|0/5|
|LowerSNR|Dev_ThinLowQ|0.876443|0.869959|0.74|0.900157|0/5|
|LowerSNR|Dev_ThinPairWideContrast|0.850296|0.840608|1.14|0.853191|0/5|
|Main|Dev_DeepContrast|0.574176|0.556713|3.04|0.604127|0/5|
|Main|Dev_ModerateQ|0.633795|0.619438|2.27|0.640172|0/5|
|Main|Dev_ThinLowQ|0.839945|0.829179|1.28|0.854938|0/5|
|Main|Dev_ThinPairWideContrast|0.785495|0.776288|1.17|0.785651|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.523442|32.95|0.523442|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.602548|-5.07|0.602548|0/1|
|StressVelocityMinus3|Dev_ThinLowQ|0.730533|0.808273|-10.64|0.808273|0/1|
|StressVelocityMinus3|Dev_ThinPairWideContrast|0.889971|0.777797|12.60|0.777797|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.523436|43.17|0.523436|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|0.602553|94.52|0.602553|0/1|
|StressVelocityPlus3|Dev_ThinLowQ|0.949055|0.808246|14.84|0.808246|0/1|
|StressVelocityPlus3|Dev_ThinPairWideContrast|0.831104|0.777797|6.41|0.777797|0/1|

未通过项目：
- StressVelocityMinus3/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressVelocityMinus3", "Model": "Dev_ModerateQ", "n": 1, "baseline": 0.573489825328423, "candidate": 0.602547689974169, "improvement_percent": -5.066849203315038, "baseline_median": 0.573489825328423, "candidate_median": 0.602547689974169, "baseline_p90": 0.573489825328423, "candidate_p90": 0.602547689974169, "baseline_p95": 0.573489825328423, "candidate_p95": 0.602547689974169, "baseline_worst": 0.573489825328423, "candidate_worst": 0.602547689974169, "baseline_failures": 0, "candidate_failures": 0}
- StressVelocityMinus3/Dev_ModerateQ/EventAmplitudeError：{"baseline": 0.456130081560017, "candidate": 0.529648338314038, "allowance": 0.02280650407800085, "n": 1}
- StressVelocityMinus3/Dev_ModerateQ/AVOInterceptError：{"baseline": 0.460863378464853, "candidate": 0.531675206649635, "allowance": 0.023043168923242652, "n": 1}
- StressVelocityMinus3/Dev_ModerateQ/AVOGradientError：{"baseline": 0.559270771693208, "candidate": 0.598050477392463, "allowance": 0.027963538584660403, "n": 1}
- StressVelocityMinus3/Dev_ThinLowQ/mean_noninferiority：{"Suite": "StressVelocityMinus3", "Model": "Dev_ThinLowQ", "n": 1, "baseline": 0.730533402330318, "candidate": 0.80827272927937, "improvement_percent": -10.641447290578698, "baseline_median": 0.730533402330318, "candidate_median": 0.80827272927937, "baseline_p90": 0.730533402330318, "candidate_p90": 0.80827272927937, "baseline_p95": 0.730533402330318, "candidate_p95": 0.80827272927937, "baseline_worst": 0.730533402330318, "candidate_worst": 0.80827272927937, "baseline_failures": 0, "candidate_failures": 0}
- StressVelocityMinus3/Dev_ThinLowQ/EventAmplitudeError：{"baseline": 0.561503630670948, "candidate": 0.771617112730198, "allowance": 0.028075181533547402, "n": 1}
- StressVelocityMinus3/Dev_ThinLowQ/AVOInterceptError：{"baseline": 0.554811128531479, "candidate": 0.766411757074745, "allowance": 0.02774055642657395, "n": 1}
- StressVelocityMinus3/Dev_ThinLowQ/AVOGradientError：{"baseline": 0.570926955950594, "candidate": 0.774780991153852, "allowance": 0.028546347797529698, "n": 1}
- main_equal_model_improvement：{"baseline": 0.7083528693454353, "candidate": 0.695404181052186, "relative_improvement": 0.018279996952951794}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
