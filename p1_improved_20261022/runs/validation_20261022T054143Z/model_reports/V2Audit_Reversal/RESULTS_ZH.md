# V2Audit_Reversal 单模型结果

仅评价已完成模型；最终分组验收和跨模型统计仍待全部完成。

| 条件 | 原版误差 | 新版误差 | 相对原版改善 | 相对 Frozen 改善 |
|---|---:|---:|---:|---:|
|Main|0.744542|0.625917|15.93%|16.97%|
|LowerSNR|0.786535|0.721643|8.25%|7.95%|
|StressNonRickerPhase|0.820291|0.703544|14.23%|15.77%|
|StressVelocityError|0.833542|0.585395|29.77%|30.12%|
|StressResidualNMO|0.811348|0.679972|16.19%|17.77%|
|StressColoredNoise|0.728846|0.629889|13.58%|16.67%|
|StressVelocityMinus3|3.528332|0.585396|83.41%|81.65%|

单模型门槛：FAIL；所有方法求解失败行数：0。负改善表示变差。

单模型未通过项目：
- StressVelocityError/V2Audit_Reversal/AVOGradientError：{"baseline": 0.7418825325840824, "candidate": 0.8944293791673911, "allowance": 0.037094126629204124, "n": 36}
