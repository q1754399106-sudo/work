# ConfidenceOnly 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 36 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 17.13%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.649913|4.85|0.720206|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.625990|13.62|0.670649|0/5|
|LowerSNR|Dev_ThinPairWideContrast|0.850296|0.724232|14.83|0.763878|0/5|
|Main|Dev_DeepContrast|0.574176|0.533097|7.15|0.615504|0/5|
|Main|Dev_ModerateQ|0.633795|0.505764|20.20|0.550803|0/5|
|Main|Dev_ThinPairWideContrast|0.785495|0.613211|21.93|0.665378|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.504541|35.37|0.504541|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.487595|14.98|0.487595|0/1|
|StressVelocityMinus3|Dev_ThinPairWideContrast|0.889971|0.710549|20.16|0.710549|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.504540|45.22|0.504540|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|0.487602|95.57|0.487602|0/1|
|StressVelocityPlus3|Dev_ThinPairWideContrast|0.831104|0.710551|14.51|0.710551|0/1|

未通过项目：
- LowerSNR/Dev_ThinPairWideContrast/WeakEventError：{"baseline": 0.792735677493171, "candidate": 0.945876250982606, "allowance": 0.03963678387465855, "n": 1}
- LowerSNR/Dev_ThinPairWideContrast/PinchoutTailError：{"baseline": 0.000696547384979592, "candidate": 0.00171764838762927, "allowance": 0.001, "n": 1}
- Main/Dev_ModerateQ/PinchoutTailError：{"baseline": 0.00196538838310368, "candidate": 0.00299182676044623, "allowance": 0.001, "n": 1}
- Main/Dev_ThinPairWideContrast/WeakEventError：{"baseline": 0.672854011465345, "candidate": 0.930714960293388, "allowance": 0.03364270057326725, "n": 1}
- Main/Dev_ThinPairWideContrast/PinchoutTailError：{"baseline": 0.000568308576057571, "candidate": 0.00173286808532089, "allowance": 0.001, "n": 1}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
