# CombinedUncertainty 对比 Frozen：效果评价

**结论：PROVISIONAL_PASS。** 配对观测 126 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 7.31%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|SoftLowQ|0.951024|0.769123|19.13|0.934132|0/12|
|LowerSNR|ThinHighFrequency|0.789360|0.743515|5.81|0.758886|0/12|
|LowerSNR|VelocityReversal|0.668010|0.616999|7.64|0.635938|0/12|
|Main|SoftLowQ|0.773073|0.675838|12.58|0.810795|0/30|
|Main|ThinHighFrequency|0.729228|0.703269|3.56|0.740695|0/30|
|Main|VelocityReversal|0.555080|0.527782|4.92|0.548386|0/30|

未通过项目：
- 无数值门槛失败；仍需完成独立性审计。

触边/目标失配诊断：{'residual target mismatch': 58}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
