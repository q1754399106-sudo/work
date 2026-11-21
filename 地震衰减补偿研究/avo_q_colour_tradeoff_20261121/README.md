# AVO 中 Q 误差与反射系数颜色的最优折中（方向 A 第 0 步）

标签：**机制研究（不立项）**。2026-11-21。

- `PLAN_ZH.md`：方案与写死的停止规则（e3df339，算数前提交）。
- `DERIVATION_ZH.md`：完整推导；第 5.3–5.4 节为事后。
- `RESULTS_STEP0_ZH.md`：S1 不过、S2 过（机制撤回）、S3 不过 → 不立项；保留的机制结论。
- `scripts/tradeoff.py`（蒙特卡罗与闭式）、`s2_stretch.py`（第二关缓存）、`s2_ideal.py`（事后理想化拉伸核对）、`verdict.py`。
- `results/`：`tradeoff.json`、`s2.json`、`s2_ideal.json`、`verdict.json`。

测井读取用本地缓存 `baselines_volve_20261029/cache/inputs/multiwell_tie_volve_20261028`（不在 git 中，可用环境变量 `E9` 改路径）。Python：codex 自带 3.12。
