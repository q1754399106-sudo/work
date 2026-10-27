# V 场景家族的新种子复验

**状态：已完成。** 结果见 [RESULTS_ZH.md](RESULTS_ZH.md)：W、X 两类都复现了翻转（标准版一侧稳定更好），支持“场景家族”解释。方案见 [PLAN_ZH.md](PLAN_ZH.md)，运行前哈希见 `PRE_RUN_HASHES.txt`。

上游：`qscan_sc_round3_20261026`（V4/V5 上标准版一侧更好）、`generator_factors_20261026`（三项可分离的生成器差别都不能复现这个翻转）。本轮用 V4/V5 的原始生成代码和全新种子，生成 6 个层状模板类场景 W1–W6 和 3 个复杂地质类场景 X1–X3，检验标准版的优势会不会在同家族的新场景上重现。

## 内容

- `scripts/snapshot_inputs.py`：round3 运行副本（`core/`、`source/`、`run_combo.m` 和只改误差档的 `run_combo3.m`）、soft_refit 的 `primary.json`、round3 的 `candidate_policy.json`、novelty 的 SC 代码 → `cache/`。
- `scripts/make_scenes.m`：生成 W1–W6、X1–X3 → `models/`、`truth/`。
- `scripts/run_matlab.py`：`make_scenes`、`stage_dev`、`run_combo3`。
- `scripts/run_sc.py`：SC、QScan＋SC。
- `scripts/score.py` → `analysis/rows.csv`、`analysis/decisions.json`。
- `scripts/plot_gathers.py` → `figures/`（运行后补画，2026-10-28）：`overview.png` 为 9 个场景的观测道集（Moderate 噪声）与真值目标；`W1.png`…`X3.png` 每个场景一张，上排 Q 正确（观测、真值、标准版、SC），下排 Q −30%（FullCalRefit、QScanRefit、QScan＋标准版、QScan＋SC）。除观测外共用真值目标的色标，标题中的相对误差逐项与 `rows.csv` 核对一致。需要 matplotlib（本地装在 `cache/plotdeps`，不上传）。

## 运行

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_matlab.py make_scenes
& $py scripts\run_matlab.py stage_dev
& $py scripts\run_matlab.py run_combo3 W1 W2 W3 W4 W5 W6 X1 X2 X3
& $py scripts\run_sc.py W1 W2 W3 ...      # 最多三个进程
& $py scripts\score.py
```

需要 MATLAB R2024b 和同级目录 `../novelty_evidence_20261025/cache/deps`（scipy，本地）。`cache/` 不上传。
