# CombinedUncertainty 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 48 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 3.13%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.650931|4.70|0.693072|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.697784|3.71|0.712232|0/5|
|LowerSNR|Dev_ThinLowQ|0.876443|0.864205|1.40|0.890933|0/5|
|LowerSNR|Dev_ThinPairWideContrast|0.850296|0.832024|2.15|0.839075|0/5|
|Main|Dev_DeepContrast|0.574176|0.542674|5.49|0.583556|0/5|
|Main|Dev_ModerateQ|0.633795|0.608067|4.06|0.622351|0/5|
|Main|Dev_ThinLowQ|0.839945|0.824348|1.86|0.845201|0/5|
|Main|Dev_ThinPairWideContrast|0.785495|0.769511|2.03|0.774361|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.770415|1.31|0.770415|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.542507|5.40|0.542507|0/1|
|StressVelocityMinus3|Dev_ThinLowQ|0.730533|0.728645|0.26|0.728645|0/1|
|StressVelocityMinus3|Dev_ThinPairWideContrast|0.889971|0.881663|0.93|0.881663|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.919984|0.11|0.919984|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|11.498683|-4.56|11.498683|1/1|
|StressVelocityPlus3|Dev_ThinLowQ|0.949055|0.946023|0.32|0.946023|0/1|
|StressVelocityPlus3|Dev_ThinPairWideContrast|0.831104|0.813057|2.17|0.813057|0/1|

未通过项目：
- StressVelocityPlus3/Dev_ModerateQ/mean_noninferiority：{"Suite": "StressVelocityPlus3", "Model": "Dev_ModerateQ", "n": 1, "baseline": 10.9972302930992, "candidate": 11.4986828986245, "improvement_percent": -4.559808171335322, "baseline_median": 10.9972302930992, "candidate_median": 11.4986828986245, "baseline_p90": 10.9972302930992, "candidate_p90": 11.4986828986245, "baseline_p95": 10.9972302930992, "candidate_p95": 11.4986828986245, "baseline_worst": 10.9972302930992, "candidate_worst": 11.4986828986245, "baseline_failures": 1, "candidate_failures": 1}
- StressVelocityPlus3/Dev_ModerateQ/AVOInterceptError：{"baseline": 2.64212577555423, "candidate": 2.80029148618854, "allowance": 0.1321062887777115, "n": 1}
- StressVelocityPlus3/Dev_ModerateQ/AVOGradientError：{"baseline": 2.33643497812194, "candidate": 2.61250536654978, "allowance": 0.116821748906097, "n": 1}
- main_equal_model_improvement：{"baseline": 0.7083528693454353, "candidate": 0.6861502583855421, "relative_improvement": 0.03134399805623689}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
