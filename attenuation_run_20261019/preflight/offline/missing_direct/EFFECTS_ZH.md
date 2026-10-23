# FixtureCandidate 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 126 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 10.00%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|SoftLowQ|0.892977|0.803679|10.00|0.978519|0/12|
|LowerSNR|ThinHighFrequency|0.796037|0.716434|10.00|0.736284|0/12|
|LowerSNR|VelocityReversal|0.675363|0.607826|10.00|0.623004|0/12|
|Main|SoftLowQ|0.763758|0.687383|10.00|1.123013|6/30|
|Main|ThinHighFrequency|0.735211|0.661690|10.00|0.687660|0/30|
|Main|VelocityReversal|0.559791|0.503812|10.00|0.542014|0/30|

未通过项目：
- matched_prior_direct_present："CSV presence only; identical prior/Q and tuning-budget provenance must be audited separately."

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
