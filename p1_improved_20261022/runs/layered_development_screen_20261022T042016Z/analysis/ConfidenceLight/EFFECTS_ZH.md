# ConfidenceLight 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 48 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 14.79%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|Dev_DeepContrast|0.683021|0.659135|3.50|0.736116|0/5|
|LowerSNR|Dev_ModerateQ|0.724652|0.639941|11.69|0.687237|0/5|
|LowerSNR|Dev_ThinLowQ|0.876443|0.814479|7.07|0.850998|0/5|
|LowerSNR|Dev_ThinPairWideContrast|0.850296|0.731657|13.95|0.785734|0/5|
|Main|Dev_DeepContrast|0.574176|0.538578|6.20|0.625047|0/5|
|Main|Dev_ModerateQ|0.633795|0.515343|18.69|0.569011|0/5|
|Main|Dev_ThinLowQ|0.839945|0.742572|11.59|0.780789|0/5|
|Main|Dev_ThinPairWideContrast|0.785495|0.617964|21.33|0.672339|0/5|
|StressVelocityMinus3|Dev_DeepContrast|0.780672|0.504541|35.37|0.504541|0/1|
|StressVelocityMinus3|Dev_ModerateQ|0.573490|0.487595|14.98|0.487595|0/1|
|StressVelocityMinus3|Dev_ThinLowQ|0.730533|0.675522|7.53|0.675522|0/1|
|StressVelocityMinus3|Dev_ThinPairWideContrast|0.889971|0.710549|20.16|0.710549|0/1|
|StressVelocityPlus3|Dev_DeepContrast|0.921016|0.504540|45.22|0.504540|0/1|
|StressVelocityPlus3|Dev_ModerateQ|10.997230|0.487602|95.57|0.487602|0/1|
|StressVelocityPlus3|Dev_ThinLowQ|0.949055|0.675513|28.82|0.675513|0/1|
|StressVelocityPlus3|Dev_ThinPairWideContrast|0.831104|0.710551|14.51|0.710551|0/1|

未通过项目：
- LowerSNR/Dev_ThinLowQ/WeakEventError：{"baseline": 1.84548601554034, "candidate": 2.415835858476, "allowance": 0.09227430077701701, "n": 1}
- LowerSNR/Dev_ThinLowQ/PinchoutTailError：{"baseline": 0.00313809667711491, "candidate": 0.00744356924381208, "allowance": 0.001, "n": 1}
- LowerSNR/Dev_ThinPairWideContrast/WeakEventError：{"baseline": 0.792735677493171, "candidate": 0.908948247360999, "allowance": 0.03963678387465855, "n": 1}
- Main/Dev_ThinLowQ/WeakEventError：{"baseline": 1.98815813969812, "candidate": 2.76351398641582, "allowance": 0.09940790698490601, "n": 1}
- Main/Dev_ThinLowQ/PinchoutTailError：{"baseline": 0.00229046924436405, "candidate": 0.00541116387815904, "allowance": 0.001, "n": 1}
- Main/Dev_ThinPairWideContrast/WeakEventError：{"baseline": 0.672854011465345, "candidate": 0.826330935724644, "allowance": 0.03364270057326725, "n": 1}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
