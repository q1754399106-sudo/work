# 更换模型后的对比结果

## Material Passport

- 类型：合成地震补偿实验；状态：实际执行完成。
- 问题：冻结优化方法在4个新地质配置中，是否仍优于最初分步法复现？
- MATLAB退出状态见 runs/changed_geology_frozen_20261023T010236Z/execution.json；原始数据、源码快照、全部结果均保留。
- 116个原有MATLAB源码文件哈希与上一轮完全一致；只增加实验入口和新模型计划。

## 设计

4个全新地质模型，每个覆盖5种场景、18/12 dB、2个新种子，合计80组普通噪声；另加20组12 dB、AR(1)=0.35有色噪声。独立校准18 dB。每个方法共用同一份观测和校准资料。
主对照FixedSplit保留原固定双岭参数1e-3；FixedSplitSharedCal额外使用与优化方法完全相同的完整校准结果。保留上一版和未剪枝结果供归因检查。
全部模型/种子/场景在运行前登记于PLAN.json，算法规则未根据本轮结果调整。

| 模型 | Q | 反射时间/s | 主频/Hz |
|---|---|---|---|
| Changed_DeepLowQ | [22, 32, 45, 75] | [0.43, 0.79, 1.04] | 50 |
| Changed_ThinPolarity | [55, 35, 90, 60] | [0.33, 0.66, 0.672] | 68 |
| Changed_LowContrast | [85, 65, 100, 130] | [0.32, 0.62, 0.96] | 38 |
| Changed_HighFreqQReverse | [38, 85, 30, 100] | [0.36, 0.7, 0.732] | 88 |

## 弱反射评分修正说明（事后复核）

薄层极性反转模型的WeakSection第6道，第二反射被场景起伏移至0.672秒，与第三反射完全重合。单道数据无法分别确定两者振幅；旧WeakEventError仍评分这一不可辨识分量，因此原始数值会误导比较。
以下表格与汇总图的弱反射误差，采用与原EventAmplitudeError相同的可辨识条件：cond(W转置W)<=1e6。对所有方法统一执行；薄层WeakSection的5组观测各排除第6道的单事件振幅评分，其余弱反射道保留。全部100组观测的完整11道波形、AVO及其他指标均未因此删除或修改。
这是运行后发现退化模板而进行的评分修正，并非运行前登记的新主指标。原评分原样保留在metrics.csv和summary.json的legacy_means；修正评分在metrics_identifiable.csv。修正未用于调参或重新筛选模型。

| 原始未修正弱反射评分 | 原分步法 | 优化方法 |
|---|---:|---:|
| 普通噪声 | 1.84284 | 1.71812 |
| 有色噪声 | 1.97283 | 1.98462 |

## 普通噪声80组

| 方法 | 波形误差 | 弱反射误差 | AVO梯度误差 | 断层邻域误差 | 尖灭尾部误差 |
|---|---:|---:|---:|---:|---:|
| FixedSplit | 0.837666 | 1.67038 | 0.45554 | 0.856729 | 0.0272893 |
| FixedSplitSharedCal | 0.809306 | 1.87769 | 0.45554 | 0.815385 | 0.026858 |
| DecoupledWeak | 0.616008 | 1.76233 | 0.522336 | 0.637454 | 0.00463256 |
| FullCalRefit | 0.202353 | 0.368818 | 0.0300209 | 0.105932 | 0.00429502 |
| PrunedFullCalRefit | 0.175495 | 0.501143 | 0.0337061 | 0.0855024 | 0 |

相对FixedSplit，优化方法波形误差降低：76/80。
未改善案例：[{"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_12dB_20262411", "baseline": 0.944170647016944, "optimized": 1.99966774762492}, {"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_12dB_20262412", "baseline": 0.931823620008453, "optimized": 1.68544427147787}, {"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_18dB_20262411", "baseline": 0.917653111700438, "optimized": 1.01829553539573}, {"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_18dB_20262412", "baseline": 0.914031387274925, "optimized": 0.934493147580974}]

相对FixedSplitSharedCal，优化方法波形误差降低：76/80。
未改善案例：[{"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_12dB_20262411", "baseline": 0.944170647016944, "optimized": 1.99966774762492}, {"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_12dB_20262412", "baseline": 0.931823620008453, "optimized": 1.68544427147787}, {"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_18dB_20262411", "baseline": 0.917653111700438, "optimized": 1.01829553539573}, {"case": "Changed_DeepLowQ_CoherentNoiseAVO_normal_18dB_20262412", "baseline": 0.914031387274925, "optimized": 0.934493147580974}]

## 有色噪声20组

| 方法 | 波形误差 | 弱反射误差 | AVO梯度误差 | 断层邻域误差 | 尖灭尾部误差 |
|---|---:|---:|---:|---:|---:|
| FixedSplit | 0.861564 | 1.81841 | 0.459509 | 0.884832 | 0.0301691 |
| FixedSplitSharedCal | 0.839426 | 2.0469 | 0.459509 | 0.851524 | 0.0303369 |
| DecoupledWeak | 0.683007 | 1.79596 | 0.603352 | 0.725901 | 0.00588559 |
| FullCalRefit | 0.292744 | 0.572504 | 0.0570713 | 0.147489 | 0.0133469 |
| PrunedFullCalRefit | 0.224811 | 0.765272 | 0.0474888 | 0.0934638 | 0 |

相对FixedSplit，优化方法波形误差降低：19/20。
未改善案例：[{"case": "Changed_DeepLowQ_CoherentNoiseAVO_colored_12dB_20262413", "baseline": 0.936170432797862, "optimized": 1.74052091208989}]

相对FixedSplitSharedCal，优化方法波形误差降低：19/20。
未改善案例：[{"case": "Changed_DeepLowQ_CoherentNoiseAVO_colored_12dB_20262413", "baseline": 0.936170432797862, "optimized": 1.74052091208989}]

## 逐模型普通噪声结果

| 模型 | 原分步波形误差 | 优化波形误差 | 原分步弱反射误差 | 优化弱反射误差 |
|---|---:|---:|---:|---:|
| Changed_DeepLowQ | 0.893219 | 0.366636 | 0.822857 | 0.111711 |
| Changed_ThinPolarity | 0.777526 | 0.140979 | 5.0062 | 1.41659 |
| Changed_LowContrast | 0.896189 | 0.087533 | 0.163955 | 0.33173 |
| Changed_HighFreqQReverse | 0.783732 | 0.10683 | 0.688511 | 0.144537 |

## 分场景普通噪声结果

| 场景 | 原分步波形误差 | 优化波形误差 |
|---|---:|---:|
| SmoothAVO | 0.752897 | 0.0547594 |
| CoherentNoiseAVO | 0.811338 | 0.524408 |
| FaultSection | 0.855396 | 0.0906001 |
| WeakSection | 0.846311 | 0.0452339 |
| PinchoutSection | 0.92239 | 0.162471 |

按“优化波形误差/原分步波形误差”自动选择的最差案例：Changed_DeepLowQ_CoherentNoiseAVO_normal_12dB_20262411，图见figures/worst_relative_case.png；不是随机或代表性样本。

## 弱反射绝对误差补充

弱反射相对误差会受真值振幅大小影响，以下补充绝对系数RMSE；不替代主指标，也不据此排除不利结果。评分在真值子波模板上估计系数，模板只用于事后评价。

| 模型 | 真弱反射系数RMS | 原分步绝对RMSE | 优化绝对RMSE |
|---|---:|---:|---:|
| Changed_DeepLowQ | 0.00659165 | 0.00542399 | 0.00073636 |
| Changed_ThinPolarity | 0.00256672 | 0.0128495 | 0.00363599 |
| Changed_LowContrast | 0.00128063 | 0.000209967 | 0.000424824 |
| Changed_HighFreqQReverse | 0.00554617 | 0.0038186 | 0.000801629 |

## 结论与退化记录

普通噪声下，弱反射平均误差高于原分步法的模型数为1/4：Changed_LowContrast。这些结果不支持“所有新模型的弱反射都改善”。
深层低Q叠加相干干扰的道集显示，约0.65秒附近非真实反射的干扰被明显放大。其真实反射位置的AVO指标仍可较小，因此需要结合完整波形和伪事件评估，不能仅以AVO指标判断结果。
本轮保持算法冻结，未根据退化案例调整阈值或重新运行选择有利结果。新增同校准分步法用于分辨校准与反演路线的作用。

## 数值与适用范围

优化阶段未收敛标志：0；最大振幅重估正规方程残差：1.49871e-14。
相对误差使用真值评分，真值和模型名称不传入反演算法。弱反射指标只在WeakSection计算，AVO指标只在两个AVO场景计算；CSV中的NaN表示不适用。尖灭振幅不可辨识处沿用原评分排除规则，整体波形仍完整评分。
本轮更换地层速度、密度、Q、反射时刻和主频，场景模板、零相位Ricker子波家族、四层/11道和正演物理保持不变。这里只是新地质配置验证，不能称为跨子波、跨正演体系或实测验证。80/20组观测来自4个地质配置，不能当作80/20个独立地质模型。
FixedSplit为按最初毕设思路实现的双岭分步法复现，不能等同于已逐行核对的毕设原始源码。
全部100幅道集共用每案例统一振幅/误差色标；预先固定展示每模型12dB第一个种子的SmoothAVO和WeakSection共8幅波形图，不按效果选图。
弱反射组筛选可能删除真实弱信号；改善百分比不等于振幅恢复比例。指标均值未附置信区间，不据此声称统计显著。

复现命令：使用scripts/run_local.py --entry run_changed_models --label replay --timeout 5400 --run-pointer <新路径>；后处理使用绘图环境运行scripts/summarize_and_plot.py。
