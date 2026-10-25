# CombinedUncertainty 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 246 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 7.31%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|DevMain|Dev_DeepContrast|0.559601|0.521686|6.78|0.567740|0/15|
|DevMain|Dev_ModerateQ|0.627602|0.593352|5.46|0.619157|0/15|
|LowerSNR|SoftLowQ|0.951024|0.769123|19.13|0.934132|0/12|
|LowerSNR|ThinHighFrequency|0.789360|0.743515|5.81|0.758886|0/12|
|LowerSNR|VelocityReversal|0.668010|0.616999|7.64|0.635938|0/12|
|Main|SoftLowQ|0.773073|0.675838|12.58|0.810795|0/30|
|Main|ThinHighFrequency|0.729228|0.703269|3.56|0.740695|0/30|
|Main|VelocityReversal|0.555080|0.527782|4.92|0.548386|0/30|
|MultiCalMain|SoftLowQ|0.841042|0.696103|17.23|0.810278|0/10|
|MultiCalMain|ThinHighFrequency|0.726272|0.700745|3.51|0.736494|0/10|
|MultiCalMain|VelocityReversal|0.550892|0.524810|4.73|0.545512|0/10|
|NewGeologyMain|Holdout_DeepLowQ|0.611691|0.588649|3.77|0.745894|0/15|
|NewGeologyMain|Holdout_HighContrastReversal|0.585150|0.554838|5.18|0.572513|0/15|
|NewGeologyMain|Holdout_StrongAtten|0.896768|0.722750|19.40|0.984859|0/15|
|NewGeologyMain|Holdout_ThinInterbed|0.629957|0.595274|5.51|0.631744|0/15|

未通过项目：
- MultiCalMain/SoftLowQ/AVOGradientError：{"baseline": 0.620129296988298, "candidate": 0.6544855183589797, "allowance": 0.0310064648494149, "n": 4}
- NewGeologyMain/Holdout_ThinInterbed/WeakEventError：{"baseline": 1.35183097947591, "candidate": 1.4390702190878033, "allowance": 0.0675915489737955, "n": 3}

触边/目标失配诊断：{'residual target mismatch': 131}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
