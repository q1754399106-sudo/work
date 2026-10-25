# CombinedUncertaintyV2 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 360 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 4.59%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|V2Audit_LowQDeep|0.864689|0.840025|2.85|1.206672|3/90|
|Main|V2Audit_LowQDeep|0.813569|0.776243|4.59|1.137883|3/90|
|StressColoredNoise|V2Audit_LowQDeep|0.855275|0.763444|10.74|0.845562|0/36|
|StressNonRickerPhase|V2Audit_LowQDeep|0.912445|0.853047|6.51|1.156279|3/36|
|StressResidualNMO|V2Audit_LowQDeep|0.853875|0.782820|8.32|1.148256|3/36|
|StressVelocityError|V2Audit_LowQDeep|15.956874|0.768510|95.18|1.161748|3/36|
|StressVelocityMinus3|V2Audit_LowQDeep|0.841736|0.768504|8.70|1.161607|3/36|

未通过项目：
- LowerSNR/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 90}
- Main/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 90}
- StressNonRickerPhase/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 36}
- StressResidualNMO/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 36}
- StressVelocityMinus3/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 36}
- main_equal_model_improvement：{"baseline": 0.8135689665312342, "candidate": 0.7762434352632185, "relative_improvement": 0.045878754971638624}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
