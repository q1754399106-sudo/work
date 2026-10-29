# 实验 14：QScan 选 Q 被压缩的机制——子波锚定与稀疏权重（2026-10-29/30）

冻结方案 `PLAN_ZH.md`，结果 `RESULTS_ZH.md`。

- 检验数据全新：构造线 10227，噪声线 10149/10273/10355；T 组 20 个（c 0.6–1.45），S 组 12 个（F-12/F-14 长测井、噪声 ×0.5/×2）。开发集为实验 12 的 20 个剖面，只用来定校正量。
- 结论：H1 子波锚定 **支持**，H2 稀疏权重决定偏差与斜率 **支持**，H3 校正后的 QScan-W 不劣于校正后的 SR **相反**，H4（压力组）无法判定。机制成立；开发集定的偏差搬不到新数据上，QScan-W 不推荐作方法。SR 不校正在新数据上最好（T 组 |对数误差| 0.046）。

| 路径 | 内容 |
|---|---|
| `scripts/run14.py` | 冻结的运行程序（dev / test 两阶段） |
| `scripts/extra14.py` | 运行中途加的辅助进程（倒序取数据集，结果相同） |
| `scripts/analyze14.py` | 冻结的校准与判定 |
| `scripts/sens14.py` | 收敛敏感性、偏差迁移、按真 c 分列（描述） |
| `scripts/fig14.py` | 汇总图 `figures/selection_T.png` |
| `scripts/gathers14.py` | 各方法选出的 c 做 SC-a 补偿后的道集图 `figures/gathers/`（描述） |
| `scripts/m1_14.py` | 道集图各面板的 M1 汇总 `analysis/gathers_m1.json`（描述） |
| `scripts/env14.py` | 导入实验 7–13 的程序（实验 13 输入快照） |
| `results/dev`、`results/test` | 每个数据集一个 JSON（选择、目标值、收敛、迭代数） |
| `analysis/` | `calibration.json`、`summary.json`、`sensitivity.json` |
| `logs/` | 运行日志 |

`cache/`（道集图的中间数组）可由 `gathers14.py` 重建，不上传。
