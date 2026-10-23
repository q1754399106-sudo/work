# ConfidenceIdentityPenalty 对比 Frozen：效果评价

**结论：PROVISIONAL_PASS。** 配对观测 24 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 12.93%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.658452|3.60|0.735932|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.640712|11.58|0.688425|0/5|
|Main|Dev_DeepContrast|0.574176|0.537658|6.36|0.623437|0/5|
|Main|Dev_ModerateQ|0.633795|0.514074|18.89|0.566730|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.504541|35.37|0.504541|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.487595|14.98|0.487595|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.504540|45.22|0.504540|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|0.487602|95.57|0.487602|0/1|

未通过项目：
- 无数值门槛失败；仍需完成独立性审计。

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
