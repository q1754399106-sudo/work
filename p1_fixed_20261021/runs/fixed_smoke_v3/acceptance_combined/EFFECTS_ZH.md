# CombinedUncertainty 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 6 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 10.13%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|Main|SoftLowQ|0.750662|0.674592|10.13|0.779944|0/5|
|StressVelocityError|SoftLowQ|7.471158|8.497215|-13.73|8.497215|1/1|

未通过项目：
- StressVelocityError/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "SoftLowQ", "n": 1, "baseline": 7.47115833239619, "candidate": 8.49721459909905, "improvement_percent": -13.733563405472339, "baseline_median": 7.47115833239619, "candidate_median": 8.49721459909905, "baseline_p90": 7.47115833239619, "candidate_p90": 8.49721459909905, "baseline_p95": 7.47115833239619, "candidate_p95": 8.49721459909905, "baseline_worst": 7.47115833239619, "candidate_worst": 8.49721459909905, "baseline_failures": 1, "candidate_failures": 1}
- StressVelocityError/SoftLowQ/EventAmplitudeError：{"baseline": 3.44756234155844, "candidate": 3.91267679304114, "allowance": 0.17237811707792203, "n": 1}
- StressVelocityError/SoftLowQ/AVOInterceptError：{"baseline": 2.65751126114753, "candidate": 3.01326155051287, "allowance": 0.1328755630573765, "n": 1}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
