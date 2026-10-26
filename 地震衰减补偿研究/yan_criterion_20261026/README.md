# Yan et al. 2023 选 Q 准则的偏移检验

- `PLAN_ZH.md`：预登记方案（提交 d203930）；`PRE_RUN_HASHES.txt`：运行前的代码、方案与输入哈希。
- `RESULTS_ZH.md`：结果与判定。
- `scripts/`：
  - `snapshot_inputs.py`：从 HEAD 复制 Git 输入（冻结的 lasso 与字典代码、G1–G8 与 H1–H6 的模型、噪声函数、归档的误差曲线与偏移汇总），并记录两份本地扫描缓存（294 个文件）的 SHA-256；
  - `yan_stats.m`：重建含噪道集（σ 须与扫描缓存逐位相同），对 13 个候选因子运行与 QScan 相同的 lasso，记录 QScan 目标 J（argmin 须等于归档的 QScan 选择）和 Yan 准则所需的逐道残差平方、Σ|r|^{1/2} → `cache/yan/`；
  - `score.py` → `analysis/selections.csv`、`analysis/offsets.csv`、`analysis/per_trace.csv`、`analysis/decisions.json`。

运行顺序（MATLAB R2024B，4 线程；Python 3.12.14）：

```
python scripts/snapshot_inputs.py
matlab -batch "addpath('scripts');yan_stats('G','1','2','3','4','5','6','7','8')"
matlab -batch "addpath('scripts');yan_stats('H','1','2','3','4','5','6')"
python scripts/score.py
```
