# CombinedUncertainty 对比 TargetMatched：效果评价

**结论：FAIL。** 配对观测 6 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 10.50%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|Main|SoftLowQ|0.753770|0.674592|10.50|0.779944|0/5|
|StressVelocityError|SoftLowQ|9.101687|8.497215|6.64|8.497215|1/1|

未通过项目：
- Main/SoftLowQ/AVOGradientError：{"baseline": 0.5887175596411816, "candidate": 0.6337073218284834, "allowance": 0.02943587798205908, "n": 2}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
