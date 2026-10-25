# 完整组合（QScan＋相干标准版）的偏差—方差实验

**状态：阶段 1、2 已完成。** 结果见 [RESULTS_ZH.md](RESULTS_ZH.md)（偏差—方差分解）与 [RESULTS_STAGE2_ZH.md](RESULTS_STAGE2_ZH.md)（选择器的改进空间）；冻结方案见 [PLAN_ZH.md](PLAN_ZH.md)、[PLAN_STAGE2_ZH.md](PLAN_STAGE2_ZH.md)。常数去偏 κ 已决定不做，见 [DECISION_KAPPA_ZH.md](DECISION_KAPPA_ZH.md)。

## 内容

- `scripts/combo_outputs.m`：在一个道集上计算 2×2 析因的四个估计量（FullCalRefit、标准版、QScanRefit、QScan＋标准版）以及 SoftGroupRefit 和 QScan＋Soft。与 `round3_20261024/run_combo.m` 逐句一致，只去掉了文件缓存。
- `scripts/cb_setup.m`：S1–S10 场景定义和种子规则（与 `qscan_bias_variance_20261024` 相同）。
- `scripts/run_stage0.m`：在归档道集上复现第三轮结果，在 r=1–2 上复现上一轮重复实验，并检查路径一致性、计时和内存。
- `scripts/run_reps.m`：阶段 1 的 worker。每个实现单独保存，已存在的会跳过，因此可以中断后续跑。
- `scripts/run_sweep.m`、`summarize_sweep.m`：阶段 2 的因子扫描与 K1–K2 判定；扫描输出在 `cache/sweep/`（不上传，约 4 MB/实现）。
- `scripts/make_figures.py`：`figures/bias_variance_2x2.png/pdf` 与 `figures/sweep_curves.png/pdf`（需要 Matplotlib，可用 `PYTHONPATH=attenuation_run_20261019/dependencies`）。
- `scripts/summarize.m`：偏差²/方差分解、区域与时窗、析因效应、交互作用、强反射振幅增益、5 批标准误，以及 J1–J5 判定。
- `scripts/diag_threads.m`：1 线程与 4 线程下相干标准版输出的差异（阶段 0 发现）。
- `scripts/hash_prerun.py`：冻结时记录方案、脚本与输入快照的哈希。
- `scripts/qscan_select.m`、`full_refit.m`、`noise_replica.m`、`estimate_sigmas.m`：与 `qscan_bias_variance_20261024/scripts` 中的文件逐字相同。
- `cache/`（不上传）：从 Git 对象逐字节取出的输入快照，以及每个实现的输出（single，约 1 MB/实现）。

## 运行

需要 MATLAB R2024b 和 Python。**内存提示**：每个 MATLAB `-nojvm` 进程实测约需 2 GB 内存。2026-10-24 晚，B0 的 6 个 worker 就是因为系统提交内存耗尽而失败的。启动前先确认 `FreeVirtualMemory` 足够。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_stage.py run_stage0
& $py scripts\run_stage.py run_reps 1 1 50     # 场景 1，实现 1–50；本轮 2 个 worker，各处理每个场景的一半
& $py scripts\run_stage.py summarize
```
