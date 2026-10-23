# 已加入：信噪比估计与自动正则化

日期：2026-10-17（Asia/Taipei）。全部结果由已有 ResearchMATLAB 会话实际运行，未启动新 MATLAB。

## 用简单的话说

现在程序会：**先看看数据有多吵，再决定补偿要多保守。**

它不再默认知道“加进去的噪声是多少”。第一步、第二步都能自动选参数；原先固定第一步参数的版本也保留，用来分清改进来自哪里。默认第二步使用与第一步一致的完整模型。

已把同一套 11 道原始含噪数据重新处理。估计信噪比在 **17.42–18.88 dB** 之间，合成设置为 18 dB：

| 炮检距（米）<br>Offset_m | 估计信噪比（分贝）<br>EstimatedSNR_dB | 第一步正则化参数<br>LambdaX | 第二步正则化参数<br>Lambda0 |
| --- | --- | --- | --- |
| 0 | 18.83 | 0.1284 | 0.001111 |
| 750 | 18.88 | 0.04899 | 0.0002547 |
| 1500 | 18.2 | 0.01175 | 9.959e-05 |

这些 λ 会随吸收算子改变；不能只看 dB 数值、用同一个经验公式给所有炮检距的 λ。两步白化后的尺度也不同，λx 与 λ0 不能直接横向比较大小。

![11 道估计信噪比和自动参数](gather/snr_parameters.png)

## 1. 噪声和信噪比怎么估计

默认从 **0.02–0.25 s** 的反射到达前时窗，利用中位数绝对偏差（MAD）估计噪声标准差：

$$\widehat\sigma=\frac{\operatorname{median}|d_\mathrm{window}-\operatorname{median}d_\mathrm{window}|}{0.67448975}.$$

这是自选噪声窗及标准高斯标定系数，不是使用干净信号挑选结果。当前第一反射位于约 0.38 s，默认窗适用于这组模型。噪声应近似平稳、高斯；用于新资料时，窗口必须确实主要包含背景噪声。

同时提供 BlindHaarMAD：不指定空白段，取相邻采样的正交差分 `(d(2k)-d(2k-1))/sqrt(2)`，再用 MAD 估计。它不需要 Wavelet Toolbox，但信号的高频成分可能污染估计，本轮把它作为对照而不是预设更优。

数据去均值后估计总功率 Pobs，噪声功率为 σhat²，信号功率取 max(Pobs−σhat²,0)：

$$\widehat{\mathrm{SNR}}=10\log_{10}\frac{\max(P_{\mathrm{obs}}-\widehat\sigma^2,0)}{\widehat\sigma^2}.$$

估计的是**观测到的含吸收信号相对于噪声**的信噪比，不是补偿后相对无吸收真值的信噪比。若无法分辨信号功率，返回 -Inf 或相应状态，不伪造正常 dB 值。

## 2. 如何据此选 λ

对每个观测和已知的正演算子 G，求

$$s_\lambda=\arg\min_s\|Gs-d\|^2+\lambda\|s\|^2,$$

再选择满足 **‖Gsλ−d‖≈σhat√n** 的 λ。这是残差原则：不要求程序把噪声也解释成地震信号。当估计信号功率为正时，这个噪声目标也可写为 √[n Pobs/(1+10^(SNRhat/10))]，所以 σhat 与 SNRhat 是同一估计链条的两种表达；无法分辨信号时保留原 σhat 和状态标记。

新增连续对数二分求根，不再只从 41 个 λ 中选最近值；初始范围为 10⁻⁷–1，必要时扩展。对相同 d、G，规定的噪声水平越大，选出的 λ 不减，这一性质已实际测试。

第一步按原始噪声选 λx。第二步考虑第一步噪声已变色，使用 BxBxᵀ 加 floor 的协方差白化，目标取 **σhat‖L⁻¹Bx‖F**，不把第一步输出重新当白噪声估计。默认第二步算子为 L⁻¹BxGt；`stage2_model='depth'` 保留原 G0 假设作对照。

因为 λx 由同一数据自适应选择，这里的协方差是冻结已选 Bx 后的 plug-in 近似，不是整个非线性自适应流程的严格无条件协方差。精确白化、一致算子及同一惩罚下，最终仍与一步反演有等价关系；这不是新颖性宣称。

## 3. 实际验证：不只展示改善的情况

主比较：3 个偏移距、3 个信噪比、3 个 Q 比例、两种正演、16 个种子；另用 8 个种子比较两步均自适应。共 **6480 条方法结果**，还有 192 条非白噪声压力试验。选参数时不使用干净信号；所有 Oracle 方法都单独标记为已知 σ 的参照。

估计精度：

| 估计方法<br>Estimator | 噪声标准差比值均值<br>SigmaRatioMean | 噪声标准差比值的标准差<br>SigmaRatioSD | 信噪比估计误差均值（分贝）<br>SNRErrorMean_dB | 信噪比估计误差标准差（分贝）<br>SNRErrorSD_dB |
| --- | --- | --- | --- | --- |
| 盲 Haar MAD 估计<br>BlindHaarMAD | 1.134 | 0.08562 | -1.06 | 0.6743 |
| 噪声窗 MAD 估计<br>WindowMAD | 0.9783 | 0.147 | 0.3288 | 1.319 |

默认 WindowMAD 本轮噪声强度均值约为真值的 0.978，信噪比平均偏差约 +0.33 dB、样本标准差约 1.32 dB。BlindHaarMAD 平均高估噪声约 13.4%，信噪比低估约 1.06 dB。不同场景复用了种子进行配对，不能把这些估计记录当作全部独立样本来证明显著性。

![噪声估计精度](figures/01_estimation_accuracy.png)

正确 Q、Column、18 dB，16 个种子的平均相对波形误差：

| 炮检距（米）<br>Offset_m | 固定参数分步补偿<br>FixedSplit | 已知噪声一步反演<br>OracleDirect | 噪声窗估计一步反演<br>EstimatedDirect | 已知噪声协方差分步（固定第一步）<br>OracleCovFixedX | 估计噪声协方差分步（固定第一步）<br>EstimatedCovFixedX | 盲估计一步反演<br>BlindDirect |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.4464 | 0.48 | 0.4707 | 0.4802 | 0.4709 | 0.508 |
| 750 | 0.7388 | 0.4988 | 0.4895 | 0.5067 | 0.4962 | 0.5406 |
| 1500 | 0.6067 | 0.5856 | 0.6023 | 0.623 | 0.6096 | 0.642 |

EstimatedDirect 表示只把噪声目标换成估计值；EstimatedCovFixedX 同样保留第一步固定 λ，只改变第二步参数。这两类对照有助于隔离“噪声估计”的影响。新 Oracle 也使用同一连续求根规则，不能把新旧网格选择的差异当作估计噪声的改善。

![信噪比、λ 与恢复误差](figures/02_snr_lambda_error.png)

两步都自动选 λ 的比较，统一使用同一批 8 个种子，避免不同样本数量混比：

| 炮检距（米）<br>Offset_m | 已知噪声一步反演<br>OracleDirect | 噪声窗估计一步反演<br>EstimatedDirect | 已知噪声双自适应（深度算子）<br>AdaptiveKnownDepth | 估计噪声双自适应（深度算子）<br>AdaptiveEstimatedDepth | 估计噪声双自适应（一致算子）<br>AdaptiveEstimatedConsistent |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.4834 | 0.47 | 0.517 | 0.5011 | 0.47 |
| 750 | 0.5033 | 0.4904 | 0.557 | 0.5381 | 0.4902 |
| 1500 | 0.588 | 0.6252 | 0.6661 | 0.6468 | 0.6176 |

AdaptiveEstimatedConsistent 是默认入口；AdaptiveEstimatedDepth 是把两步参数自适应后仍用 G0 的对照。后者偏差较大，说明自动选参数不能取代算子一致性。

把全部主实验等权汇总，并按 Q 比例分组：

| Q 模型比例<br>QScale | 已知噪声一步反演<br>OracleDirect | 噪声窗估计一步反演<br>EstimatedDirect | 已知噪声协方差分步（固定第一步）<br>OracleCovFixedX | 估计噪声协方差分步（固定第一步）<br>EstimatedCovFixedX | 盲估计一步反演<br>BlindDirect |
| --- | --- | --- | --- | --- | --- |
| 0.8 | 0.8278 | 1.035 | 0.786 | 0.8283 | 0.7712 |
| 1 | 0.5249 | 0.5291 | 0.5398 | 0.5289 | 0.566 |
| 1.2 | 0.6607 | 0.6517 | 0.6675 | 0.659 | 0.6864 |

这个表包含噪声估计导致的退化条件。特别是 Q 不准确、噪声被低估时，残差目标可能过小，λ 随之过小，出现过补偿。本轮全部主网格平均 EstimatedDirect 约 **0.7387**，OracleDirect 约 **0.6711**；EstimatedCovFixedX 约 **0.6720**，OracleCovFixedX 约 **0.6644**。这些总体平均采用自选等权场景，不代表实测资料的分布，也不能用某一组更好的结果宣称全面优越。

![Q 偏差下的反例](figures/03_q_mismatch.png)

## 4. 对非白噪声的适用边界

64 个种子分别测试白高斯噪声、AR(0.7) 相关噪声，以及高斯背景加 3 个幅度为 12σ 的随机脉冲：

| 随机噪声<br>Noise | 估计 / 真实噪声标准差<br>SigmaRatio | 白噪声检查通过比例<br>WhiteAssumptionOK | 相对波形误差<br>RelativeError |
| --- | --- | --- | --- |
| 有色噪声 AR(0.7)<br>ColoredAR07 | 0.9936 | 0 | 0.5574 |
| 高斯噪声＋脉冲<br>GaussianPlusImpulses | 0.7796 | 0.9688 | 0.6856 |
| 高斯白噪声<br>WhiteGaussian | 1.016 | 0.9688 | 0.5056 |

WindowMAD 对相关噪声的边际噪声强度估计接近 1，但白化仍没有使用真实相关协方差，不能认为统计模型已正确。噪声窗 lag-1 检查在本轮检出了相关噪声，但它不是完整的高斯性或平稳性检验。

脉冲噪声下 MAD 主要估计背景噪声，低估了包含脉冲在内的总噪声约 22%；lag-1 检查通常无法发现这个问题。这些情形需要进一步增加相关噪声协方差或稳健残差损失，不能仅依靠一个 SNR 数字解决。

![非白噪声压力试验](figures/04_noise_stress.png)

## 5. 已执行的检查与代码入口

- MATLAB 主验证 PASS，6480 条方法结果；行/列向量输入一致，全零输入保持全零。
- 将配置中的 `snr_db` 改成 -Inf、seed 改成 NaN，给定同一观测的输出完全不变：默认算法没有读取合成真值。
- 将数据整体乘 7，估计 σ 同比变化、SNR 不变、输出也乘 7，相对误差约 **6.90e-14**。
- 连续选择器通过矩形算子与正规方程解对照，以及规定噪声增大时 λ 单调不减检查。
- 实际跑通共享 Engine 的 `--action snr` 入口；首次集成发现官方 Engine 的 stdout 只接受 io.StringIO，已修复并重新运行通过。
- 用户配置保持，代码更改前快照在 `work/backups/snr_update_20261017/`；本次没有清空 base 工作区、改变 MATLAB 许可或新开实例。

生产接口：`snr_compensate(d, offset, c)`；默认估计噪声并选择两步参数，不需要输入真实 SNR。

可选：`opts.method='haar'` 使用自动 Haar 对照；`opts.stage1_mode='fixed'` 保留第一步参数；`opts.stage2_model='depth'` 采用历史第二步；`noise_window` 可按实际噪声段设置。`sigma_override` 仅为明确标记 Oracle 的验证接口，默认不设置。

`run_rebuild_snr` 已实际处理 11 道观测，结果在 `gather/`。`validate_snr_pipeline` 为全量实验驱动；`estimate_trace_snr` 与 `discrepancy_select` 可单独复用。旧 `run_rebuild`、原参数文件及上轮结果保持作为基础参照。

## 理论来源

MAD 的高斯标定依据 [MathWorks wnoisest 官方说明](https://www.mathworks.com/help/wavelet/ref/wnoisest.html)，本实现自行计算中位数，不调用额外工具箱。残差原则及 Tikhonov 参数选择可参见 [DTU 的 Michael Jacobsen 博士论文中的参数选择讨论](https://www2.imm.dtu.dk/pubdb/edoc/imm3380.pdf)。这些是已有原则，本次工作是将数据驱动的噪声估计接入当前吸收补偿，并如实验证其局限。
