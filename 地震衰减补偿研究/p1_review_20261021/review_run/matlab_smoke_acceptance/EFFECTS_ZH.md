# CombinedUncertainty 对比 Frozen：效果评价

**结论：PROVISIONAL_PASS。** 配对观测 5 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 10.18%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|Main|SoftLowQ|0.750662|0.674260|10.18|0.779819|0/5|

未通过项目：
- 无数值门槛失败；仍需完成独立性审计。

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
