# ConfidenceOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 24 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 14.00%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.649913|4.85|0.720206|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.625990|13.62|0.670649|0/5|
|Main|Dev_DeepContrast|0.574176|0.533097|7.15|0.615504|0/5|
|Main|Dev_ModerateQ|0.633795|0.505764|20.20|0.550803|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.504541|35.37|0.504541|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.487595|14.98|0.487595|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.504540|45.22|0.504540|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|0.487602|95.57|0.487602|0/1|

未通过项目：
- Main/Dev_ModerateQ/PinchoutTailError：{"baseline": 0.00196538838310368, "candidate": 0.00299182676044623, "allowance": 0.001, "n": 1}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
