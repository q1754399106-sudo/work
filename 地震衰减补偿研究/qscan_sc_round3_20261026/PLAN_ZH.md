# QScan＋SC 在一维数据 F3/F4、V4/V5 上的正交对照（预登记）

登记时间：2026-10-26。在任何本轮 SC 求解或评分之前写入。

## 问题

`qscan_sc_20261026` 在二维有限差分 H1–H6 上判定：Q 低估 30% 时，QScan＋SC 稳定优于 QScan＋相干标准版（6/6）。最终报告第 9 节注明它“还没有在其他生成器上确认”。

用户 2026-10-26 提供的外部验证方案把这个对照列为第一阶段。本轮把同一组合放到 round3 的一维数据上，按 {名义 Q, QScan} × {SC, 标准版, FullCalRefit} 和 7 档 Q 误差（−50% 到 +50%）做正交对照，回答：

1. **Q 有误时**，QScan＋SC 是否仍稳定优于 QScan＋标准版？分两种生成器分别判定：
   - **Kennett 组**（F3、F4）：一维 Kennett 细层正演。SC 的参数在同一生成器的 F1/F2 上选定，SC 有生成器优势。
   - **V 组**（V4、V5）：`gen_layered_scene`、`gen_complex_model`，即相干规则开发时所用的生成器。标准版有生成器优势。
2. **Q 正确时**，SC 上叠加 QScan 是否变差（与 `qscan_sc_20261026` 的 H2 同一问题）。

V 组是对 QScan＋SC 最不利的数据：Q 正确时，标准版在 V4/V5 上 4/4 胜过 SC（`novelty_evidence_20261025`）。

## 登记前已经知道的内容

- round3 的 56 行（F3/F4/V4/V5 × Moderate/Strong × 7 档）上，标准版、QScan＋标准版、FullCalRefit、QScanRefit 的输出和误差都已看到过（`round3_20261024/analysis/combo_metrics.csv`）。
- Q 正确（Q+0）时 SC 在这 8 行上的误差已看到过（`novelty_evidence_20261025/analysis/comparison.csv`）：F3/F4 上 SC 4/4 胜标准版，V4/V5 上标准版 4/4 胜。
- F3/F4、V5 都已多次复用，不是未见数据（`final_report_20261025/ERRATA_ZH.md` 第 10 条）。
- **从未计算过**：Q 有误时的 SC（名义 Q），以及任何一行的 QScan＋SC。

## 数据与方法

- **数据**：`round3_20261024/dev/{F3,F4,V4,V5}.mat` 中的 Moderate、Strong 两个含噪道集（每行一个噪声实现），评分目标取 `*_truth.mat` 的 `target`。F1/F2 不用：SC 在上面调过参。
- **Q 误差档**：e ∈ {−50, −30, −10, 0, +10, +30, +50}%，所用 Q = Q0 × (1 + e/100)，与 round3 的 `run_combo.m` 相同。
- **QScan 的选择**：取 `round3_20261024/results/<道集>_Q<e>_<噪声>.mat` 中归档的 `scan.factor` 和 `scanned.Q`，也就是 QScan＋标准版用的同一个因子。不重新运行 QScan。
- **SC**：`novelty_evidence_20261025` 提交 ff95a69 的 `benchmark.solve`、`benchmark.predictor`、`common.dictionary`，逐字使用：
  - 参数冻结为 a = 0.125、mu 相对强度 = 10（`analysis/selection.json`）；
  - KKT ≤ 0.02，最多 5000 次迭代；
  - 字典用道集的 `config`，只把 Q 换成所用 Q（SC）或 `scanned.Q`（QScan＋SC）。
- **已归档的对照**（从 round3 结果文件读取，不重算）：标准版、QScan＋标准版、FullCalRefit、QScanRefit。
- **一致性检查**：Q+0 时本轮重算的 SC 应与 `novelty_evidence_20261025/cache/solved/<道集>_<噪声>_0.125_10.npz` 一致（相对差 ≤ 1e-8），否则停止并记录。
- 环境：Python 3.12.14、numpy 2.3.5、scipy 1.18.1，与 novelty、sc_neutral、qscan_sc 相同；BLAS 单线程；最多 3 个进程。

共 112 次 SC 求解：8 个道集 × 7 档 × {名义 Q, QScan}。

## 指标

- **主指标**：相对波形误差 ‖X − T‖/‖T‖，全部道、全时段，与 novelty、round3 的 `overall` 相同。一行一个值。
- **次指标**：归一化平方误差。
- **用时（只作描述）**：SC 每次求解的秒数（Python，单线程）；round3 日志中 MATLAB 每行的用时（4 线程，包含 QScan 扫描、两次重拟合和相干收缩）。两边的语言、线程数和所含步骤都不同，不能据此判定速度优劣。

## 判定

“A 稳定优于 B”沿用 sc_neutral、qscan_sc 的三条门槛，同时满足才成立：
1. 平均误差 A ≤ 0.95 × B；
2. 至少 75% 的行（向上取整）中 A 更好；
3. 任一行 A/B ≤ 1.05。

每个判定都在两个方向上各算一次，都不成立记为“无稳定优势”。

- **R1-K（Q 有误，Kennett 组）**：F3/F4 × 2 种噪声 × 6 个非零误差档 = 24 行，QScan＋SC 对 QScan＋标准版（需 ≥ 18/24）。
- **R1-V（Q 有误，V 组）**：V4/V5 同样 24 行。
- **R1 结论规则**：
  - R1-K、R1-V 都判“QScan＋SC 稳定更好”→ **在一维生成器上确认**，与 H1–H6 一致；
  - 只有 R1-K 成立 → **依赖生成器**：在标准版的开发生成器上没有确认；
  - 任一组反方向成立 → **与 H1–H6 结论冲突**，在该生成器上 QScan＋标准版更好；
  - 其余 → 不确认。
- **R2（Q 正确）**：Q+0 的 8 行，SC 对 QScan＋SC（需 ≥ 6/8）。结论措辞与 `qscan_sc_20261026` 的 H2 相同。

## 只描述

- 分组、分误差档、分道集的平均与胜负；
- QScan＋SC 对 SC（名义 Q）的挽回幅度；QScan＋SC 对 QScanRefit；
- 用户 2026-10-26 提供的外部验证方案（阶段 1）中的两条通过标准，按原文核对：“QScan＋标准版平均误差优于或持平 QScan＋SC”；“运行耗时比 SC 降低 50% 以上”（后者只能给出带上述限制的描述）；
- 次指标；收敛情况。

## 收敛与失败处理

与 sc_neutral、qscan_sc 相同：未收敛的输出照常计入判定并标出，另报只用收敛行的敏感性分析；不加大迭代次数、不改容差。程序出错需要重跑时，记录原因，不改设计。

## 不做

- 不重新运行或修改 QScan，不为任何方法调参，不用 F1/F2、U1–U4。
- 看到任何本轮结果后，不改分组、门槛或指标。

## 边界

- 已见数据：除 QScan＋SC 和 Q 有误时的 SC 外，其余数字登记前都已看到。
- 统计单位少：每组只有 2 个模型、4 个道集，同一道集的 6 个误差档高度相关；每行一个噪声实现。24 行不能当作独立样本，不报告显著性检验。
- 两组数据对两种方法都不中立（见“问题”），所以分组判定，不合并成一个总判定。
- SC 是 FISTA 独立实现，不是作者原版代码。QScan 的选择基于 FullCalRefit 的 lasso 目标，不是针对 SC 设计的。
