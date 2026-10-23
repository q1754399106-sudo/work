# CalShrink050Light 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 48 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 11.18%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.661815|3.10|0.729793|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.664111|8.35|0.711760|0/5|
|LowerSNR|Dev_ThinLowQ|0.876443|0.830726|5.22|0.878630|0/5|
|LowerSNR|Dev_ThinPairWideContrast|0.850296|0.760608|10.55|0.832100|0/5|
|Main|Dev_DeepContrast|0.574176|0.544375|5.19|0.614696|0/5|
|Main|Dev_ModerateQ|0.633795|0.544935|14.02|0.600112|0/5|
|Main|Dev_ThinLowQ|0.839945|0.765709|8.84|0.821703|0/5|
|Main|Dev_ThinPairWideContrast|0.785495|0.661552|15.78|0.741584|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.504541|35.37|0.504541|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.487595|14.98|0.487595|0/1|
|StressVelocityMinus3|Dev_ThinLowQ|0.730533|0.675522|7.53|0.675522|0/1|
|StressVelocityMinus3|Dev_ThinPairWideContrast|0.889971|0.710549|20.16|0.710549|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.504540|45.22|0.504540|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|0.487602|95.57|0.487602|0/1|
|StressVelocityPlus3|Dev_ThinLowQ|0.949055|0.675513|28.82|0.675513|0/1|
|StressVelocityPlus3|Dev_ThinPairWideContrast|0.831104|0.710551|14.51|0.710551|0/1|

未通过项目：
- LowerSNR/Dev_ThinLowQ/WeakEventError：{"baseline": 1.84548601554034, "candidate": 2.00961626987388, "allowance": 0.09227430077701701, "n": 1}
- Main/Dev_ThinLowQ/WeakEventError：{"baseline": 1.98815813969812, "candidate": 2.19775339700736, "allowance": 0.09940790698490601, "n": 1}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
