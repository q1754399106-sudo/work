# CombinedUncertainty V2：冻结候选使用说明

是否通过最终确认，应以同一运行目录的 `analysis/RESULTS_ZH.md`、`analysis/acceptance.json` 及完整性审计为准。开发筛选通过不能替代完整确认。

## 调用冻结版本

本轮冻结运行目录：

`C:/Users/qq175/Documents/ChatGPT/研究/p1_improved_20261022/runs/validation_20261022T054143Z`

```matlab
addpath('C:/Users/qq175/Documents/ChatGPT/研究/p1_improved_20261022/runs/validation_20261022T054143Z/source/core_p1');
out = combined_uncertainty_v2(D, offsets, public, calibrationD, calibrationOffsets, 'Section');
S = out.gather;                 % 当前默认交付 Split 补偿波形
S_direct = out.direct;          % 共享校准与结构先验配置的 Direct 参考
S_split = out.split;
```

上例假定五个输入变量已由实际观测及公开背景模型准备好。输入是 MATLAB 数值数组和结构体，不是文件名。完整实验期间不应修改冻结目录中的方法、策略或阈值。

## 输入与已验证范围

| 输入 | 内容 |
|---|---|
| `D` | 时间样本×11道的观测矩阵；当前确认配置为601样本、采样间隔0.002秒 |
| `offsets` | 11个已知偏移距，单位米，必须属于 `public.offset` 网格 |
| `calibrationD` | 同一背景地层模型的独立校准观测，时间样本×11道；不是无噪声真值 |
| `calibrationOffsets` | 校准观测对应的11个已知偏移距 |
| `public` | 仅包含允许的采集、背景模型及固定算法配置字段，见下文 |
| 域 | `'Offset'` 表示偏移距道集，`'Section'` 表示横向剖面；由采集方式指定，不由真实测试场景标签选择 |

本轮验证使用 `public.offset = 0:150:1500`。Offset观测采用这11个偏移距；Section观测为11道、每道偏移距均为750米，校准观测仍为完整的11偏移距道集。其他采样、道数、采集布置和真实资料条件不在本轮确认范围内。

`public` 允许字段为：

```text
dt, nt, t, offset, vp, h, fref, nfft,
lambda_x, lambda_fixed, lambda_grid, cov_floor, nominal_frequency
```

时间单位为秒，速度为米/秒，厚度和偏移距为米，频率为Hz。`vp` 有四层，`h` 给出前三层厚度。背景界面深度由厚度累加得到，预期零偏移到时由 `cumsum(2*h./vp(1:3))` 得到，因此背景界面位置属于已知信息。层数、速度和厚度不是由本方法从零恢复。

`nominal_frequency` 是事先给定的名义主频，实际主频从校准观测估计。初始子波形式采用Ricker；压力试验检验相位等偏离，但这不等于已验证任意未知子波。

真实Q、真实主频、横波速度、密度、真实事件到时、反射系数和无噪声波形不得放入 `public`；额外字段会被生产入口拒绝。若从合成配置生成公开结构体，应使用白名单函数 `uq_public_config`，而不是直接传入包含真值的完整配置。

## 候选改变了什么

先用校准观测拟合有效Q、子波主频及有限范围内的整体速度修正，额外自由度由未参与训练的校准道检验；再尝试四层Q校准，同样通过校准道交叉验证决定是否接受。该过程不读取补偿评分或真实Q。

冻结策略 `CalShrink025Light` 在Offset域使用完整的新校准。在Section域，对有效Q校准与分层校准的点估计取几何收缩：

\[
Q_i^*=Q_{\mathrm{fixed}}^{0.75}Q_{\mathrm{layer},i}^{0.25},\qquad
f^*=f_{\mathrm{fixed}}^{0.75}f_{\mathrm{layer}}^{0.25}.
\]

整体速度修正保留。Section结构约束采用噪声置信度下限0.2、模态权重混合0.25、残差目标倍数1、结构权重3。每个域的参数均在开发阶段确定，Direct与Split共享校准参数、结构算子H、噪声估计及选参尺度；Q误差模态权重则按各自传播算子的奇异基构造，因此不应将其称为完全相同物理正则矩阵的对照。输出仍经过完整衰减补偿求解，没有原始观测直通。

收缩是旨在控制点估计方差的确定性规则，不是经验证的贝叶斯模型后验平均。Q和频率的既有不确定性信息用于正则化；分层Q区间覆盖率、各层速度/厚度的联合反演及外部实测泛化尚未由本轮证明。

## 输出与诊断

| 字段 | 含义 |
|---|---|
| `gather`、`direct`、`split` | 补偿波形及两条路线 |
| `estimated_Q` | 本次域实际采用的Q/频率估计，可能包含 `layer_Q`；以字段内容为准 |
| `base_calibration_estimate` | Section收缩前的新校准点估计 |
| `estimated_velocity_scale` | 实际采用的整体速度修正系数，不是逐层速度反演结果 |
| `calibration_diagnostics` | 校准训练/验证划分、接受状态、拟合及边界诊断 |
| `direct_record`、`split_record`、`records` | 最终解及选参过程的求解记录 |
| `policy` | 从冻结JSON读取的完整候选策略与开发来源 |

`compensation_applied=true` 仅表示执行了补偿，不表示效果已通过验收。是否合格须结合全覆盖、全部方法的求解质量、局部误差、交叉bootstrap区间、源码哈希和实际重放结果判断。

## 证据解释

开发阶段为四模型、单校准和单观测种子，候选Main等模型平均误差较原Combined改善约6.77%，其中一个模型小幅退化约0.56%。候选通过了预先保留的工程门槛；其余失败配置与诊断实验均保留。

完整确认包含旧九模型回归与四个提前封存的新合成模型、三校准、六观测种子及二十场景/条件。旧九模型不是盲测；新四模型仍使用同一合成模拟器，也不能表述为外部实测验证。完整确认的实际结论由结果报告给出。

图形及辅助软件的确切依赖、工作流程说明与文献见 `source/scripts/figure_requirements.lock.txt` 和 `source/scripts/FIGURE_SOFTWARE_NOTES_ZH.md`。
