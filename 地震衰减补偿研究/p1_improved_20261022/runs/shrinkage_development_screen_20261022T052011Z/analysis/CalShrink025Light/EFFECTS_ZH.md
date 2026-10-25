# CalShrink025Light 对比 Frozen：效果评价

**结论：PROVISIONAL_PASS。** 配对观测 48 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 9.70%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.663117|2.91|0.726250|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.674914|6.86|0.720851|0/5|
|LowerSNR|Dev_ThinLowQ|0.876443|0.838105|4.37|0.888846|0/5|
|LowerSNR|Dev_ThinPairWideContrast|0.850296|0.769992|9.44|0.845134|0/5|
|Main|Dev_DeepContrast|0.574176|0.545691|4.96|0.607341|0/5|
|Main|Dev_ModerateQ|0.633795|0.560004|11.64|0.611125|0/5|
|Main|Dev_ThinLowQ|0.839945|0.776705|7.53|0.840298|0/5|
|Main|Dev_ThinPairWideContrast|0.785495|0.676264|13.91|0.766289|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.504541|35.37|0.504541|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.487595|14.98|0.487595|0/1|
|StressVelocityMinus3|Dev_ThinLowQ|0.730533|0.675522|7.53|0.675522|0/1|
|StressVelocityMinus3|Dev_ThinPairWideContrast|0.889971|0.710549|20.16|0.710549|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.504540|45.22|0.504540|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|0.487602|95.57|0.487602|0/1|
|StressVelocityPlus3|Dev_ThinLowQ|0.949055|0.675513|28.82|0.675513|0/1|
|StressVelocityPlus3|Dev_ThinPairWideContrast|0.831104|0.710551|14.51|0.710551|0/1|

未通过项目：
- 无数值门槛失败；仍需完成独立性审计。

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
