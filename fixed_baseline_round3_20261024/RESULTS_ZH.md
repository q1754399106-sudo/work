# 最初分步法在 F3/F4、V5 上的补跑

## Material Passport

- 日期：2026-10-24；MATLAB R2024b 实际运行（`matlab -batch`，51 s）。
- 目的：总结中缺少“最初方法”和最新方法在冻结后新数据上的直接对比。此前最初分步法只在 20 个开发场景上有结果。
- 方法：最早归档的 FixedSplit。使用 `round3_20261024/source/fixed_v1_{operators,apply}.m`，与 `soft_refit_20261023/source` 中的版本仅换行符不同。调用方式同 `q_dense_gathers_20261023/run_dense.m`：`fixed_v1_operators(used, repmat(offset,1,nx))`，offset 与 round3 的 `run_combo.m` 相同（`solver_offset`，默认 450）。
- 数据：round3 的 `dev/F3、F4、V5.mat`，观测、σ 和公开配置不变，Q 误差网格 −50…+50% 与 round3 相同。真值只在 `score_fixed.m` 中读取。
- 检查：
  - 正规方程残差 ≤1.2e-14，重放差 ≤5e-15（见 `analysis/fixed_checks_F3_F4_V5.json`）。
  - 同一评分函数重算 round3 的 FullCalRefit，与 `guard_metrics_*.csv` 的最大相对差为 7.4e-15，评分口径一致（见 `analysis/score_checks_F3_F4_V5.json`）。
- 运行前哈希：`logs/PRE_RUN_HASHES.txt`。
- 性质：这是补充对照，不是预登记实验，不设门槛。

## 结果

整体相对误差，越低越好。“最坏”指每个场景、噪声组合在 7 档 Q 误差中的最大值，再取平均。

| 数据 | 方法 | 平均误差 | Q 正确时 | 最坏 | 优于最初方法 |
|---|---|---:|---:|---:|---:|
| F3/F4（28 场景） | 最初分步法 | 1.051 | 0.884 | 1.730 | — |
| | FullCalRefit | 0.875 | 0.548 | 2.331 | 22/28 |
| | 相干标准版 | 0.713 | **0.458** | 1.789 | 26/28 |
| | QScanRefit | 0.536 | 0.546 | 0.555 | 28/28 |
| | **QScan + 相干标准版** | **0.468** | 0.472 | **0.479** | **28/28** |
| V5（14 场景） | 最初分步法 | 0.918 | 0.704 | 1.640 | — |
| | FullCalRefit | 1.331 | 0.605 | 4.727 | 10/14 |
| | 相干标准版 | 1.276 | 0.597 | 4.463 | 11/14 |
| | QScanRefit | 0.564 | 0.563 | 0.589 | 14/14 |
| | **QScan + 相干标准版** | **0.540** | **0.543** | **0.564** | **14/14** |

## 解读

- **QScan + 相干标准版相对最初方法：**
  - F3/F4：平均误差降 55%（1.051 → 0.468），Q 正确时降 47%，最坏情况降 72%；
  - V5：平均误差降 41%，Q 为标定值时降 23%；
  - 42 个场景全部更好。
- **不加 QScan 时有反例**：Q 低估 30–50% 时，FullCalRefit 和相干标准版会失控，误差比最初方法更大。
  - V5 上它们的平均误差（1.33 / 1.28）反而高于最初方法（0.918）；
  - F3/F4 在 Q −50% 时也是最初方法更好。
  - 在 Q 正确时，新方法始终更好。
- 所以，相对最初方法的稳定收益需要 QScan 这一步：它让新方法在 Q 有误时也不崩溃。

## 文件

- `run_fixed.m` → `results/*.mat`：每个场景的输出道集，LFS 存储。
- `score_fixed.m` → `analysis/fixed_metrics_F3_F4_V5.csv`、`score_checks_*.json`。
- 运行：在 MATLAB 中执行 `run_fixed; score_fixed`。需要同级目录 `round3_20261024`，稀疏检出时需先恢复。
- 日志：`logs/run_fixed.log`。
