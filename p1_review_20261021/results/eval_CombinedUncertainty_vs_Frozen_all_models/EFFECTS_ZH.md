# CombinedUncertainty 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 246 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 6.07%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|DevMain|Dev_DeepContrast|0.555567|0.519675|6.46|0.569001|0/15|
|DevMain|Dev_ModerateQ|0.623227|0.592813|4.88|0.616389|0/15|
|LowerSNR|SoftLowQ|0.950922|0.845410|11.10|1.225678|2/12|
|LowerSNR|ThinHighFrequency|0.789360|0.753187|4.58|0.769010|0/12|
|LowerSNR|VelocityReversal|0.668010|0.629730|5.73|0.649918|0/12|
|Main|SoftLowQ|0.773074|0.695747|10.00|0.940063|0/30|
|Main|ThinHighFrequency|0.729228|0.705976|3.19|0.741977|0/30|
|Main|VelocityReversal|0.555080|0.530855|4.36|0.548914|0/30|
|MultiCalMain|SoftLowQ|0.895024|0.750888|16.10|1.045061|2/10|
|MultiCalMain|ThinHighFrequency|0.736741|0.714065|3.08|0.756193|0/10|
|MultiCalMain|VelocityReversal|0.550469|0.527592|4.16|0.545030|0/10|
|NewGeologyMain|Holdout_DeepLowQ|0.608064|0.587504|3.38|0.747855|0/15|
|NewGeologyMain|Holdout_HighContrastReversal|0.585202|0.557535|4.73|0.573154|0/15|
|NewGeologyMain|Holdout_StrongAtten|0.884106|0.753961|14.72|1.209277|3/15|
|NewGeologyMain|Holdout_ThinInterbed|0.624575|0.592093|5.20|0.631847|0/15|

未通过项目：
- NewGeologyMain/Holdout_ThinInterbed/WeakEventError：{"baseline": 1.4179925427020568, "candidate": 1.5124526792886868, "allowance": 0.07089962713510284, "n": 3}

触边/目标失配诊断：{'residual target mismatch': 142}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
