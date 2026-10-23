# 实验方法 / Experiment methods

研究问题是误差预算怎样改变有效Q补偿中的稳定性与保幅。研究对象是带限波形，不是逐层Q恢复。所有模型采用四层速度/密度/Q、三个反射界面、601采样点、2ms采样、11道和0:150:1500m炮检距。保留原始精确分层射线与时变算子，不能假定完整非稳态算子可精确分解为两个算子的乘积。Shuey二项、无倍波/透射/几何扩散的模型仍是主要适用边界。

The target is the band-limited unattenuated waveform. Four-layer, three-event models retain the original exact layered-ray, time-varying forward operators. The full nonstationary operator is used explicitly; the experiment never assumes an exact product decomposition. These are chosen synthetic layer parameters, not independently collected field geology.

开发模型与四个测试模型在 `experiment_spec.json` 中完整写出。开发/测试ID不重叠，且入口逐项拒绝与归档历史模型相同的vp/Q/event_t三元组。新模型包括交替层Q、深部强吸收、薄低速层与阻抗抵消。这个设计保证参数组不重复；它不能单独证明真实地质泛化，也不能把同一模型中的道、场景或噪声次数当作独立地质样本。

Main18 and LowerSNR12 each contain SmoothAVO, FaultSection, PinchoutSection, and WeakSection with white observation noise. CoherentNoiseAVO is retained in its separate StressCoherent suite. This pre-execution deviation from the original five-scene per-suite checklist avoids counting an exact duplicate of SmoothAVO when coherent contamination is absent; all five scene types remain in the full experiment. Each held-out model has three separately synthesized calibration noise realizations and six independent observation realizations nested within each calibration. The same nested seeds are reused across methods, scenes, and SNRs for paired comparisons, never across calibration and observation namespaces.

校准和反演只接收D、炮检距与公开运动学参数白名单。真Q、真实主频、真实振幅与干净目标只在 `p1_synthesize` 和后验评分中使用。错速度是公开速度整体偏5%，厚度与真实射线保持原合成模型；残余NMO是逐事件符号交替、远炮检距最大12ms的二次时移，无NMO拉伸。源相位为30°，形状压力为延迟9ms且频率0.68倍、权重0.32的Ricker双波，不属于单一相位Ricker模板。Q压力是逐层±25%的交替Q扰动。噪声压力分别为AR(1)系数0.65、干净道集范数25%的相干波、以及0.4%采样点中额外10倍名义噪声标准差的脉冲。除所列单一因素外均保持同一模型、种子与SNR。

All methods share the same observation, public input, calibration data, and point Q/frequency estimate. Baseline and Frozen reuse their original algorithm. Candidate branches separate uncertainty in the noise scale, model error, and their combination. ScalarBudgetDamping, NoiseCovarianceOnly, and NoStructureCombined are predeclared mechanism controls. Each branch produces Direct and Split with the identical estimated prior/covariance information and selected scale. PriorHash and SharedInformationHash allow route consistency checks. Normal residuals, all solver flags, discrepancy endpoints, and missed targets are retained.

开发不搜索超参数。候选内部基于观测的差异原则搜索预算在配置中固定；校准bootstrap固定8次，噪声bootstrap固定64次。冻结封存参数、模型划分、种子、超参数预算、估计器、工程阈值及源码。测试不能回馈调参。

Acceptance follows the provisional engineering goals in NEXT_EXPERIMENT: at least5% equal-model Main waveform improvement over Frozen Split, no more than2% suite/model mean degradation, bounded event/weak/AVO errors, bounded pinchout leakage, no larger error>1 failure fraction, and normal-equation relative residual at most2e-7. The same-prior Direct results are reported. Hierarchical descriptive resampling treats models and then calibrations as sampling levels; the four-model interval does not establish significance. All means, medians, p90/p95, worst cases, paired differences, failures, and adverse local outcomes must be reported. None of these alone authorizes replacing Frozen.

预定示例固定第一测试模型/第一校准/第一观测，固定道1/6/11，不能按结果挑选好看的道集。`worst_cases/` 按每方法/路线最大有限RelativeError的固定规则在线保存JSON/MAT，包含有限求解器失败、真值、观测及同观测Frozen Direct/Split。最多16个案例，非有限失败保留CSV和日志。入口不在MATLAB中画图，评价阶段在MATLAB外根据这些真实保存波形绘制最差案例。这个后验规则只服务报告，绝不进入选参。
