# V2Audit_LowQDeep 单模型结果

仅评价已完成模型；最终分组验收和跨模型统计仍待全部完成。

| 条件 | 原版误差 | 新版误差 | 相对原版改善 | 相对 Frozen 改善 |
|---|---:|---:|---:|---:|
|Main|0.788382|0.776243|1.54%|4.59%|
|LowerSNR|0.846179|0.840025|0.73%|2.85%|
|StressNonRickerPhase|0.877060|0.853047|2.74%|6.51%|
|StressVelocityError|18.262642|0.768510|95.79%|95.18%|
|StressResidualNMO|0.809640|0.782820|3.31%|8.32%|
|StressColoredNoise|0.823588|0.763444|7.30%|10.74%|
|StressVelocityMinus3|0.788597|0.768504|2.55%|8.70%|

单模型门槛：FAIL；所有方法求解失败行数：0。负改善表示变差。

单模型未通过项目：
- LowerSNR/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 90}
- Main/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 90}
- StressNonRickerPhase/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 36}
- StressResidualNMO/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 36}
- StressVelocityMinus3/V2Audit_LowQDeep/failure_count：{"baseline": 0, "candidate": 3, "n": 36}
- main_equal_model_improvement：{"baseline": 0.8135689665312342, "candidate": 0.7762434352632185, "relative_improvement": 0.045878754971638624}
