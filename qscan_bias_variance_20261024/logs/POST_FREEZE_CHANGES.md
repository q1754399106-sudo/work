# 冻结后的改动记录

`logs/PRE_RUN_HASHES.json` 记录于 2026-10-24T04:20:14Z，即阶段 1（无噪声对照）和阶段 2（重复噪声）开始之前。之后的改动如下，均未改变求解器、噪声生成、评分口径或预先判定：

1. `scripts/summarize.m`：派生的逐样本均值/方差矩阵改存到本地 `cache/`（可由 `results/rep_*.mat` 精确重建），不再写入 `results/`。数值计算未变。改动发生在阶段 1 已完成、阶段 2 运行中，summarize 尚未运行。
2. `scripts/analyze.py`：在看到阶段 1 无噪声结果之后、阶段 2 重复噪声结果产生之前编写；D1–D5 的判定式按 `PLAN_ZH.md` 原文实现。
3. `scripts/summarize.m`：增加可选场景子集参数，用于在 F1/Moderate 仍在运行时先检查汇总流程（`summarize(2,3,4)`）。默认仍汇总全部四个场景，计算未变。最终表格由全部场景的运行生成。
4. `summarize_snr` 首次运行：14:14 已打印 `DONE snr summary` 并写出 `analysis/snr_probe.csv`，但 MATLAB 进程未退出，14:25 由我终止，因此该次执行记录为 FAILED（`logs/summarize_snr_attempt1*`）。随后在负载较低时重跑。
