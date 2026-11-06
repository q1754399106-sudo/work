# 式 17 为什么没有被锚定：判据、求解器还是数据

接 `qscan_orig_solver_20261106`。方案 `PLAN_ZH.md`（运行前提交）。

- `scripts/repo_e17.py`：旧流程（与 `qscan_anchor_theory_20261030/scripts/crit_kappa.py` 相同的解）上加算式 17。
- `scripts/scan_q2.m`、`run_irls.m`：Wang 2021 IRLS 上的式 17、chengD、C2。
- `scripts/wang_op.py`、`lasso_scan.py`、`run_lasso.py`：仓库 lasso 求解器放到 Wang 算子上。
- `scripts/gen_cache.m`：生成两个求解器共用的数据（`cache/`，可重算，不上传）。
- `scripts/analyze.py`：评分，写入 `analysis/factor.json`。
- `scripts/dev_*`：开发用（种子 901–902），不进入判定。

依赖：运行前要把 `qscan_orig_solver_20261106/scripts` 和 `qscan_anchor_theory_20261030/scripts` 恢复到工作区。
