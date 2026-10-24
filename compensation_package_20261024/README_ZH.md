# 吸收补偿（Q 补偿）推荐版本：统一入口

本目录把截至 2026-10-24 的推荐方案整理成一个 MATLAB 函数 `q_compensate`，代码与已归档实验逐值一致（见“自检”）。

## 怎么调用

```matlab
addpath('compensation_package_20261024');
out = q_compensate(D, cm);                    % 相干标准版，Q 按给定值
out = q_compensate(D, cm, 'mode', 'qscan');   % 先用 QScan 选整体 Q 因子，再做相干标准版
```

- `D`：nt×nx 共偏移距剖面道集；全零道视为坏道。
- `cm`：标定配置，需含 `t, dt, nt, nfft, fref, fpeak, vp(1×4), h(1×3), Q(1×4)`；偏移距默认 `cm.solver_offset`，否则 450 m，可用 `'offset'` 指定。
- `'sigma'`：每道噪声尺度；默认用 0.02–0.25 s 窗口 MAD 估计，该窗口有信号时请自行提供。
- `'rule'`：`'standard'`（默认，λ 0.5/16）、`'v1'`（λ 1/8，弱反射更保守）、`'soft'`（不用邻道信息）。

输出 `out.gather` 为补偿结果；另含未收缩结果 `unshrunk`、原软收缩 `soft`、实际使用的 Q、QScan 因子、逐组相干度 c、λ、增益、坏道回退标记和数值检查。

## 选哪种模式

| 情况 | 推荐 | 依据 |
|---|---|---|
| Q 来自可信的独立参考（VSP、测井） | `'standard'` | 独立正演 F3/F4 上 Q 正确时误差为 FullCalRefit 的 0.84 倍 |
| Q 只来自地面资料标定，可能有误 | `'qscan'` | Q 误差 ±50% 时最坏误差约为 FullCalRefit 的 0.12–0.20 倍；复杂模型标定 Q 下比标准版再低 9% |
| 需要最大限度保留弱反射 | `'rule','v1'` | 标准版弱反射误差比 v1 约多 1.6% |

已知限制：

- QScan 在细层（非稀疏反射系数）模型上偏向少补偿约 10–28%，原因是 lasso 稀疏先验与密集反射系数不匹配；改用其他选择准则或分层因子都没有更好（`qscan_criterion_20261024`、`round4_20261024`、`qcriteria_20261024`）。因此 Q 已经正确时，`'qscan'` 约差 3%。
- 线性相干干扰与真反射一样连续，相干收缩不会压制它。
- 复杂模型的主要剩余误差来自子波 / 模型失配，收缩规则改善有限（约 1%）。
- 所有证据来自合成数据（同源卷积正演与独立 Kennett 细层正演），尚无实测资料验证；单道、零相位 Ricker 字典、已知层位几何的假设仍然成立。

## 自检

```matlab
run('compensation_package_20261024/tests/run_selftest.m')
```

用随包的 F3、F4、V4 观测，经本入口重算后与 `round3_20261024/results` 中归档的“标准版”和“QScan + 标准版”道集逐值比较（相对差 0），QScan 因子完全一致；另检查全坏道输入输出为零、λ=1 时与原软收缩一致。

## 来源

- 相干软收缩：`coherence_soft_20261023`（v1，验收通过）、`coherence_soft_v2_20261024`（标准版 λ 0.5/16）。
- QScan：`bae_q_20261024`；组合实验：`round3_20261024`。
- `source/` 为冻结的求解器文件，与上述实验相同；`core/coherence_scores.m`、`core/apply_gain.m` 来自 v2 / round3。
- 示例：`examples/demo.m`（F4 强噪声、Q 低估 30%）。
