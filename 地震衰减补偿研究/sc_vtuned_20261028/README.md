# SC 在 V 家族上调参后的公平复比

**状态：已完成。** 结论“部分解释”，见 [RESULTS_ZH.md](RESULTS_ZH.md)。方案见 [PLAN_ZH.md](PLAN_ZH.md)（预登记提交 43ab07f），运行中的修改见 [RUN_NOTES_ZH.md](RUN_NOTES_ZH.md)。

上游：`vfamily_20261027`。它的边界一节指出无法区分两种解释：相干规则适配了 V 家族，或者 V 家族本身有利于标准版。原因是标准版的规则在 V1–V3 上开发，而 SC 的参数是在 Kennett F1/F2 上选的。本轮在 V1–V3 上给 SC 重新选参（25 组候选，标准版阶段 A 为 12 组），再在 W1–W6、X1–X3 上用同样的门槛复比。

## 内容

- `scripts/run.py`：`dev` 在 V1–V3 上跑候选网格（读真值，属开发数据）；`eval` 在 W/X 上跑 SC_tuned 和 QScan＋SC_tuned（不读真值）。求解器原样调用 `novelty_evidence_20261025/scripts/benchmark.py`，并校验其哈希。
- `scripts/choose.py`：汇总网格、执行一次边界扩展、冻结 `analysis/selection.json`。
- `scripts/score.py`：W/X 判定。标准版一侧用 `vfamily_20261027/analysis/rows.csv` 的归档误差。
- `analysis/`：`selection.json`、`extension.json`、`rows.csv`、`decisions.json`。
- `logs/`：求解日志。

## 运行

纯 Python（numpy、scipy），不需要 MATLAB。需要 LFS 文件：`coherence_soft_20261023/inputs/validation/V1–V3.mat` 与 `truth/`，以及 `vfamily_20261027/models/`、`truth/`。

```
python scripts/run.py dev
python scripts/choose.py
python scripts/run.py dev --extended
python scripts/choose.py --final
python scripts/run.py eval
python scripts/score.py
```

云端 2 核耗时：开发约 4.5 小时（其中约 2 小时花在后来出局、未收敛的候选上），评估约 1 小时。`cache/` 不上传。
