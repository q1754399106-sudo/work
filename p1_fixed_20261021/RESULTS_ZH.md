# MATLAB 修改后重跑结果

进程状态：COMPLETED，总耗时 262.3 秒；共 108 行指标。

全部记录的最大正规方程残差：9.998e-08；残差或标志失败记录：0。

执行开始后源码哈希变化：0 个文件。

| 场景 | Frozen 误差 | 修改后 Combined 误差 | 改善 |
|---|---:|---:|---:|
| CoherentNoiseAVO | 1.132113 | 0.779944 | +31.11% |
| FaultSection | 0.666569 | 0.670918 | -0.65% |
| PinchoutSection | 0.625761 | 0.612496 | +2.12% |
| SmoothAVO | 0.667976 | 0.641961 | +3.89% |
| WeakSection | 0.660892 | 0.667639 | -1.02% |

五场景平均相对改善 10.13%；单模型、单校准、单观测种子，不能据此证明普遍改善。

| 方法 | 修改前五场景均值 | 修改后五场景均值 |
|---|---:|---:|
| Baseline | 0.876420 | 0.876420 |
| CombinedFixedTarget | — | 0.692213 |
| CombinedUncertainty | 0.674260 | 0.674592 |
| CombinedUncertainty_NoStructure | 0.792904 | 0.791706 |
| Frozen | 0.750662 | 0.750662 |
| NoiseUncertaintyOnly | 0.673960 | 0.673960 |
| QModelErrorOnly | 0.746660 | 0.797934 |
| StrongerDamping | 0.768211 | 0.768211 |
| TargetMatched | — | 0.753770 |

## 修正后的压力测试与验收

实际速度扰动 +3%，厚度保持不变；重建后算子差异 Frobenius 范数 2.337315。

Frozen 相对误差 7.471158，Combined 相对误差 8.497215，候选较 Frozen 退化 13.73%。两者误差都大于 1。

该条件下 Q=15.014928，Q 触边=True，校准优化器 exitflag=0。求解器收敛不等于校准可靠或波形恢复成功。

这是速度与相应垂直走时/校准失配的联合压力条件，不是保持全部事件走时不变的纯偏移距动校正实验。

acceptance_combined：**FAIL**，14/17 项通过。
- StressVelocityError/SoftLowQ/mean_noninferiority：{"Suite": "StressVelocityError", "Model": "SoftLowQ", "n": 1, "baseline": 7.47115833239619, "candidate": 8.49721459909905, "improvement_percent": -13.733563405472339, "baseline_median": 7.47115833239619, "candidate_median": 8.49721459909905, "baseline_p90": 7.47115833239619, "candidate_p90": 8.49721459909905, "baseline_p95": 7.47115833239619, "candidate_p95": 8.49721459909905, "baseline_worst": 7.47115833239619, "candidate_worst": 8.49721459909905, "baseline_failures": 1, "candidate_failures": 1}
- StressVelocityError/SoftLowQ/EventAmplitudeError：{"baseline": 3.44756234155844, "candidate": 3.91267679304114, "allowance": 0.17237811707792203, "n": 1}
- StressVelocityError/SoftLowQ/AVOInterceptError：{"baseline": 2.65751126114753, "candidate": 3.01326155051287, "allowance": 0.1328755630573765, "n": 1}

acceptance_target_matched：**FAIL**，16/17 项通过。
- Main/SoftLowQ/AVOGradientError：{"baseline": 0.5887175596411816, "candidate": 0.6337073218284834, "allowance": 0.02943587798205908, "n": 2}

回归检查全部通过：True。

新增目标匹配对照用于检验收益来源，不能只凭 StrongerDamping 对照排除所有阻尼解释。

方法为不确定性平均的经验模态正则化；未宣称实现完整数据误差协方差。校准协方差包含经验下限。

完整九模型、多校准、多种子网格未在本次执行；历史扩展集的两个局部失败不能因本次小范围回归通过而宣布消失。

修订版采用独立 MATLAB 实现；未对原包缺失依赖的 Python 数值双生补造验证。
