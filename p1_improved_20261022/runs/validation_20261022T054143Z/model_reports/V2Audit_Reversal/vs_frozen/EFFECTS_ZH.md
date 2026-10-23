# CombinedUncertaintyV2 对比 Frozen：效果评价

**结论：FAIL。** 配对观测 360 组，路线 Split。
主试验等模型平均波形误差相对变化：改善 16.97%（负数表示退化）。

这是工程门槛检查；尚需核对新地质独立性、校准重复、信息公平性与统计区间。通过不自动替换主版。

| Suite / 组别 | Model / 模型 | Baseline / 原版 | Candidate / 候选 | Improvement / 改善% | Worst / 候选最差 | Failures / 候选误差>1 |
|---|---|---:|---:|---:|---:|---:|
|LowerSNR|V2Audit_Reversal|0.783931|0.721643|7.95|0.771138|0/90|
|Main|V2Audit_Reversal|0.753878|0.625917|16.97|0.671710|0/90|
|StressColoredNoise|V2Audit_Reversal|0.755936|0.629889|16.67|0.648973|0/36|
|StressNonRickerPhase|V2Audit_Reversal|0.835298|0.703544|15.77|0.721500|0/36|
|StressResidualNMO|V2Audit_Reversal|0.826899|0.679972|17.77|0.696966|0/36|
|StressVelocityError|V2Audit_Reversal|0.837721|0.585395|30.12|0.610906|0/36|
|StressVelocityMinus3|V2Audit_Reversal|3.190497|0.585396|81.65|0.610909|0/36|

未通过项目：
- StressVelocityError/V2Audit_Reversal/AVOGradientError：{"baseline": 0.7418825325840824, "candidate": 0.8944293791673911, "allowance": 0.037094126629204124, "n": 36}

触边/目标失配诊断：{}。
完整检查见 acceptance.json；所有逐观测配对差保存在 paired_differences.csv。
