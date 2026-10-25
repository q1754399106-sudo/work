# Q 不确定性感知反射簇实验结果

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-10-23
- Verification Status: ANALYZED; actual MATLAB execution and numerical limiting controls passed; no independent reproduction
- Version Label: q_cluster_v1_completed

## 结论

首轮预登记工程验收：**未通过**。新模型等权正常道误差相对 FullCalRefit 改善 -2.40%，低能量区误差改善 +0.91%。负数表示退化。
这是新机制的实际首轮检验，不是通过重命名 FullCalRefit 得到的方法。即使某一项指标改善，也不能跳过整体、局部与消融门槛。未证明首创性。

## 方法与创新假设

保持原独立校准和 L1 支持发现，使用校准 logQ 的局部误差尺度构造固定先导模型误差白化；对每个反射簇分解组合模式，再按固定复杂度目标降低簇秩并联合重估振幅。目标是保留可辨识的组合波形，而非整簇删除。
创新候选集中在 Q 误差驱动的局部自由度选择。SVD、固定先导误差协方差、自适应参数化和稀疏反演本身均不是本次首创。完整规则见 PLAN_ZH.md，全部候选源码、开发版本快照和冻结哈希保留。

## 实际执行范围

开发为原复杂模型×2条件×22道=44条波形；验证为3个新地质配置×2条件×22道=132条波形；Q点估计扰动另44条波形。共220条观测/扰动处理条件，每条比较6种恢复方法，共1320条方法波形；重复扰动不是新观测或独立样本。数值控制另计。
三个新配置分别改变层位密度、Q深度趋势、主频、相位、偏移距、地质随机种子、校准噪声与观测噪声。依然继承同族准1D属性驱动正演，不是完整二维粘弹性正演或实测验证。
预先从99道中均匀选22道。本次选中位置恰好不含原生成器的第23/78坏道，因此全道和正常道评分相同，不能从本轮推断坏道鲁棒性。

## 冻结验证均值（每模型、每条件等权）

| 方法 | 波形误差 | 低能量区误差 | 深部误差 | 相关系数 |
|---|---:|---:|---:|---:|
| FullCalRefit | 0.558942 | 0.907118 | 0.632591 | 0.827884 |
| PrunedFullCalRefit | 0.574538 | 0.981918 | 0.654353 | 0.814184 |
| GlobalTSVD | 0.561024 | 0.914758 | 0.635220 | 0.825922 |
| ClusterModesNoQ | 0.570582 | 0.904824 | 0.646826 | 0.818686 |
| QAwareRidge | 0.555461 | 0.911416 | 0.629026 | 0.830248 |
| QAwareCluster | 0.572338 | 0.898836 | 0.648189 | 0.818763 |

## 逐模型核对与校准诊断

| 模型 | FullCalRefit误差 | 候选误差 | 候选相对改善 | Q点估计触边 |
|---|---:|---:|---:|---|
| FreshFoldLowQ | 0.673980 | 0.693335 | -2.87% | False |
| FreshDenseQReverse | 0.520095 | 0.532698 | -2.42% | True |
| FreshThinLowFrequency | 0.482751 | 0.490982 | -1.71% | False |

其中 Q 触边记录仍纳入全部评分；没有按校准或处理效果删除不利模型。触边提示校准受限，不等于候选已自动修正偏差。

## 配对消融

- QAwareCluster 相对 PrunedFullCalRefit 的波形误差改善：+0.38%。
- QAwareCluster 相对 GlobalTSVD 的波形误差改善：-2.02%。
- QAwareCluster 相对 ClusterModesNoQ 的波形误差改善：-0.31%。
- QAwareCluster 相对 QAwareRidge 的波形误差改善：-3.04%。

Q误差机制是否有收益，应同时查看 QAwareRidge 对 FullCalRefit，以及 QAwareCluster 对 ClusterModesNoQ；簇秩机制应查看 QAwareCluster 对 QAwareRidge。不能把组合结果都归功于 Q。

## 预登记门槛

| 门槛 | 结果 |
|---|---|
| overall_improves_2pct | FAIL |
| each_model_within_2pct | FAIL |
| better_than_noQ_1pct | FAIL |
| better_than_Qridge_1pct | FAIL |
| low_energy_within_2pct | PASS |
| numerical | PASS |
| all_pass | FAIL |

候选在6个模型×条件组合中的整体误差胜出 0/6；条件来自3个模型，不把该胜率解释成独立成功概率。

## Q点估计扰动

第一新模型的同一 Moderate 观测，名义 Q 与0.8/1.25倍估计 Q；使用相同的名义校准协方差来源，不使用真实 Q。以下最差误差与输出变化分别衡量准确性和敏感性，稳定不等于正确。

| 方法 | 三点最差波形误差 | 最大输出变化 |
|---|---:|---:|
| FullCalRefit | 0.832157 | 1.254496 |
| PrunedFullCalRefit | 0.833389 | 1.291846 |
| GlobalTSVD | 0.833435 | 1.265911 |
| ClusterModesNoQ | 0.835161 | 1.292421 |
| QAwareRidge | 0.818282 | 1.101596 |
| QAwareCluster | 0.832512 | 1.219561 |

## 数值、流程与解释边界

- 开发、验证和敏感性 MATLAB 运行均实际完成；验证/敏感性最大正规方程相对残差 6.392e-15。全部数值检查详见各运行 checks.json。
- 校准协方差为带固定弱先验和尺度截断的局部线性化估计，不是已验证的贝叶斯后验/置信区间；第四层可能缺乏校准信息。
- 使用单个固定先导计算 Q 模型误差，忽略更高阶项、源子波空间变化和有色噪声协方差；不能称为完整模型误差似然。
- 模式方向、支持和误差权重依赖观测先导，复杂度项不是经校准的显著性检验。更小拟合目标并不保证真值误差更低。
- 普通截断、NoQ消融和候选共用已有 L1 支持。候选无法恢复支持发现阶段已遗漏且未补入的反射。
- validation_metrics.csv 中 ObservedReplayResidual 字段实际为 ||输出-输入||/||输入||，仅表示处理变化量，不能解释为正演重放残差或去噪改善。本报告不将其用于晋级。
- 模式路径中的数据拟合项可由 objective_path 减去 sigma² log(nt)×rank 得到；白化不同的方法其数值不能直接当作相同似然比较。
- 首轮开发已显示整体退化，仍按已登记规则完成验证以检验取舍；未根据新模型结果调参。开发快照之后仅增加缺字段防护、数值控制与预登记Q扰动入口，未改变已有开发条件下估计器。
- 本次冻结验证的三个地质配置现已被查看。若利用本轮结果修改算法，这些配置应转为开发/回归资料，不能继续充当后续独立确认集。
- 11/11类统计谬误检查覆盖见 analysis/validation_RESULTS_ZH.md；不做显著性、独立样本置信区间或跨地区泛化声明。

## 后续定位（尚未实验）

当前模式选择目标控制的是加权观测拟合与自由度，没有直接保护补偿输出 W r 的振幅。应先在开发资料上量化被删模式对输出波形的贡献，再研究带输出保真约束的模式选择。该解释是下一步假设；本轮不能证明某一个阈值或物理因素是唯一失败原因。

## 可复现入口

依赖 MATLAB R2024b（未使用额外工具箱）；数据生成/评分/绘图使用 Python 3.12、NumPy、Matplotlib、Pillow。各运行 SOURCE_MANIFEST.json 保存准确源码快照，FROZEN_CANDIDATE.json 保存验证前源码、输入和真值哈希。
重新执行已冻结输入时，运行 python run_experiment.py development 或 python run_experiment.py validation，会创建新的运行目录并保留历史。运行后 python analyze.py development / validation 生成对应评分。新一轮重复实验不要覆盖历史真值或再次运行有防覆盖保护的生成器。
先运行 prepare_sensitivity.py，再运行 run_experiment.py sensitivity 和 summarize_sensitivity.py 可复核已登记 Q 扰动；本次输入已存在，重复时直接使用现有输入运行 sensitivity，不再准备覆盖。最后运行 build_report.py 生成本报告和图件。

## 文件

- RESULTS_ZH.md：本报告。
- analysis/validation_metrics.csv：逐模型、逐条件、逐方法指标。
- analysis/validation_trace_metrics.csv：全部逐道误差。
- analysis/validation_modes.csv：支持数、保留秩、复杂度目标和数值残差。
- analysis/sensitivity_metrics.csv：所有 Q 扰动指标。
- figures/：全部PNG/PDF、来源与显示尺度清单。
- runs/：原始MATLAB结果、校准诊断、日志、PID、退出码和源码快照。