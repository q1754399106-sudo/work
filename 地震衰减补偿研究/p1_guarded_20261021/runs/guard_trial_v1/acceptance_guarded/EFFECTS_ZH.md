# GuardedNoise 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 14 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 9.92%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|Main|SoftLowQ|0.743973|0.670170|9.92|0.767339|0/10|
|StressVelocityMinus3|SoftLowQ|0.750653|1.000895|-33.34|1.000973|2/2|
|StressVelocityPlus3|SoftLowQ|8.007510|1.000895|87.50|1.000973|2/2|

未通过项目：
- StressVelocityMinus3/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityMinus3", "Model": "SoftLowQ", "n": 2, "baseline": 0.75065288041706, "candidate": 1.0008946401813752, "improvement_percent": -33.33654826253139, "baseline_median": 0.75065288041706, "candidate_median": 1.0008946401813752, "baseline_p90": 0.7557543590925536, "candidate_p90": 1.000957057503651, "baseline_p95": 0.7563920439269903, "candidate_p95": 1.0009648596689356, "baseline_worst": 0.757029728761427, "candidate_worst": 1.00097266183422, "baseline_failures": 0, "candidate_failures": 2}
- StressVelocityMinus3/SoftLowQ/failure_count：{"baseline": 0, "candidate": 2, "n": 2}
- StressVelocityMinus3/SoftLowQ/EventAmplitudeError：{"baseline": 0.709522710358759, "candidate": 0.997668500515875, "allowance": 0.035476135517937955, "n": 2}
- StressVelocityMinus3/SoftLowQ/AVOInterceptError：{"baseline": 0.6872197787080555, "candidate": 0.991792876406848, "allowance": 0.03436098893540278, "n": 2}
- StressVelocityMinus3/SoftLowQ/AVOGradientError：{"baseline": 0.719473474394835, "candidate": 0.997829253072292, "allowance": 0.03597367371974175, "n": 2}
- candidate_compensation_coverage：{"compensated_observations": 10, "total": 14, "note": "Pass-through rejection is not successful attenuation compensation."}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
