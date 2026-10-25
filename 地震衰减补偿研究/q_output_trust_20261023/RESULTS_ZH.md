# 第二轮：输出保幅信赖域实验结果

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-10-23
- Verification Status: ANALYZED; actual MATLAB runs and numerical controls completed; no independent reproduction
- Version Label: q_output_trust_v2

## 结论

预登记验收：**未通过**。正常道波形误差相对FullCalRefit改善 +0.18%，低能量区误差改善 +0.32%。负值表示退化。候选在6个模型×条件组合胜出 3/6；这些条件只来自3个地质配置。
开发按既定规则从0.5/1/2三个kappa中选定 **0.5**，新验证与Q扰动全程固定。实际发生模式删除的验证道占 100.0%；平均受保护输出方向数 59.30；相对完整Q加权参考的输出变化均值 0.0890。

## 机制与先前诊断

先前4个已查看模型的8个条件中，簇秩删除相对QAwareRidge均使真值均方误差增加。用delta=S_cluster-S_ref核对误差变化恒等式 2<error_ref,delta>+||delta||²，完整逐条件结果见analysis/previous_round_output_diagnosis.csv；仅作开发诊断。
本轮保持Q误差白化与簇内模式表示，但每次删除还需满足输出空间的累计信赖域：在参考投影超过2倍近似标准差的方向，整个候选相对完整参考的变化不能超过kappa倍该标准差。协方差通过支持到补偿波形的映射传播，约束相对完整参考，防止多步漂移。数学定义见PLAN_ZH.md和core/qtrust_refit.m。
参考可能有偏；支持选择、噪声假设、校准协方差、固定先导和尺度截断均为近似，因此不是经覆盖率验证的置信区间，也不保证真实保幅。自适应参数化、SVD和误差传播各自已有研究；本实现的首创性尚未确认。

## 冻结验证均值

| 方法 | 波形误差 | 低能量区误差 | 深部误差 | 相关系数 |
|---|---:|---:|---:|---:|
| FullCalRefit | 0.576492 | 0.925051 | 0.665522 | 0.808805 |
| PrunedFullCalRefit | 0.600056 | 1.001868 | 0.699095 | 0.784564 |
| QAwareRidge | 0.575235 | 0.928694 | 0.664983 | 0.809026 |
| QAwareCluster | 0.584580 | 0.918717 | 0.674975 | 0.802633 |
| OutputTrustNoQ | 0.577316 | 0.917792 | 0.666626 | 0.808283 |
| QAwareOutputTrust | 0.575447 | 0.922093 | 0.665198 | 0.809131 |

## 逐模型与配对消融

| 模型 | FullCalRefit误差 | 候选误差 | 相对改善 | 校准Q触边 |
|---|---:|---:|---:|---|
| V2DenseLowQ | 0.723050 | 0.723233 | -0.03% | False |
| V2BroadContrast | 0.568089 | 0.562739 | +0.94% | False |
| V2ThinQReverse | 0.438336 | 0.440370 | -0.46% | False |

本轮全部模型均纳入评分。

- 相对 QAwareCluster 的波形误差改善：+1.56%。
- 相对 QAwareRidge 的波形误差改善：-0.04%。
- 相对 OutputTrustNoQ 的波形误差改善：+0.32%。

## 预登记门槛

| 门槛 | 结果 |
|---|---|
| overall_2pct | FAIL |
| each_model_within_2pct | PASS |
| better_than_Qridge_1pct | FAIL |
| better_than_NoQ_05pct | FAIL |
| low_energy_within_2pct | PASS |
| nontrivial_deletions_10pct | PASS |
| numerical | PASS |
| all_pass | FAIL |

## Q点估计扰动

第一新模型的同一Moderate观测，Q点估计乘0.8/1/1.25；保持名义校准协方差来源，不使用真Q。

| 方法 | 三点最差误差 | 最大输出变化 |
|---|---:|---:|
| FullCalRefit | 0.799815 | 1.068226 |
| PrunedFullCalRefit | 0.810390 | 1.151738 |
| QAwareRidge | 0.794883 | 1.013164 |
| QAwareCluster | 0.808323 | 1.033547 |
| OutputTrustNoQ | 0.800244 | 1.076857 |
| QAwareOutputTrust | 0.796722 | 1.016498 |

## 补充条件协方差核验

固定支持、固定先导、已知高斯噪声与模型误差的独立稠密计算，512次模拟方差相对解析方差的偏差中位数为 3.45%，最大为 8.19%，通过原定数值要求。该检查验证条件线性误差传播，不验证真实Q误差分布或选择后的置信覆盖率。
补充检验首次失败：稠密协方差约9e-17的非对称舍入误差使重根特征分解进入一般非对称路径，得到非正交特征基；第二次增加诊断输出时暴露微小复数导致JSON写入失败。随后仅在补充检验中显式对称化协方差，保留原样本数和判据，通过核验。主估计器使用薄SVD白化，未改动。失败日志与空诊断文件均保留，详见analysis/covariance*。

## 范围、核查与后续边界

- 开发44条波形条件×7方法；新模型验证132条×6方法；额外Q扰动44条×6方法。共220条观测/扰动条件，1364条方法波形；其中真实不同的输入观测只有176条。
- 三个新地质配置与上一轮均不同，改变层密度、Q趋势、主频、相位、偏移距以及校准/观测种子；仍为同族准1D属性驱动正演，不能称为跨正演或实测验证。
- 固定22道空间抽样不含原第23/78坏道；不据此推断坏道鲁棒性。低能量区分位定义只在真值评分程序中使用，不解释为单个弱反射系数误差。
- 初始支持及残差补入沿用FullCalRefit；本轮不解决已遗漏事件、任意子波、未知层位或不可靠Q校准。
- 所有基线共用观测和校准；真值从求解输入中隔离。原MATLAB依赖哈希保持一致，开发/验证/扰动均保存实际执行源码和输入快照。
- 在独立MATLAB进程中从准备缓存重放一个选定开发道，输出相对差为0；该检查只核对缓存和估计器一致性。
- 在开发快照后只新增了预登记Q扰动入口；名义开发估计器qtrust_refit.m未改变。验证前完整冻结，未根据新模型结果改参数。
- 完整记录正规方程残差、保护方向数、累计预算使用量、删除秩、被保护规则拒绝的提案和所有开发候选。数值通过与科学效果晋级分开判断。
- 11/11类统计谬误按本描述性范围核查：分模型展示避免聚合掩盖；不从区域均值推断单事件；披露生成模型选择与先导依赖；不按结果删案例或控制后验筛选变量；基率不适用；不把开发变化当重复测量收益；保留失败；多候选仅在开发选择、不作显著性；完整冻结分析分支；归因限配对机制；受控计算不涉及观察性反向因果。
- 若继续利用本轮结果改算法，本轮三个新配置应转为开发/回归资料，后续需要新的确认资料。

## 复现与文件

MATLAB R2024b；Python 3.12 + NumPy/Matplotlib/Pillow。run_round2.py development/validation/sensitivity 各自新建运行目录；analyze_round2.py development/validation 评分；prepare_sensitivity.py（输入已有时跳过）和summarize_sensitivity.py处理扰动；build_report.py生成报告图件；audit_completion.py核查交付。
RESULTS_ZH.md为本报告，SELECTED.json保留开发选择，FROZEN_CANDIDATE.json保留验证前源码/输入/真值哈希。analysis/含全部逐道和逐条件CSV。runs/保留日志、退出码、精确快照、原始输出以及不含真值的prepared.mat缓存。图件统一尺度、不逐道归一化、不截幅。