# 完整组合（QScan＋相干标准版）的偏差—方差实验

**状态：设计稿与脚本已写好，阶段 0 尚未运行。** 设计见 [DESIGN_ZH.md](DESIGN_ZH.md)。阶段 0 通过后将改名为 `PLAN_ZH.md` 并冻结。

## 内容

- `scripts/combo_outputs.m`：在一个道集上计算 2×2 析因的四个估计量（FullCalRefit、标准版、QScanRefit、QScan＋标准版）以及 SoftGroupRefit 和 QScan＋Soft。与 `round3_20261024/run_combo.m` 逐句一致，只去掉了文件缓存。
- `scripts/cb_setup.m`：S1–S10 场景定义和种子规则（与 `qscan_bias_variance_20261024` 相同）。
- `scripts/run_stage0.m`：在归档道集上复现第三轮结果，在 r=1–2 上复现上一轮重复实验，并检查路径一致性、计时和内存。
- `scripts/run_reps.m`：阶段 1 的 worker。每个实现单独保存，已存在的会跳过，因此可以中断后续跑。
- `scripts/summarize.m`：偏差²/方差分解、区域与时窗、析因效应、交互作用、强反射振幅增益、5 批标准误，以及 J1–J5 判定。
- `scripts/qscan_select.m`、`full_refit.m`、`noise_replica.m`、`estimate_sigmas.m`：与 `qscan_bias_variance_20261024/scripts` 中的文件逐字相同。
- `cache/`（不上传）：从 Git 对象逐字节取出的输入快照，以及每个实现的输出（single，约 1 MB/实现）。

## 运行

需要 MATLAB R2024b 和 Python。**内存提示**：每个 MATLAB `-nojvm` 进程大约需要 1 GB 提交内存。2026-10-24 晚，B0 的 6 个 worker 就是因为系统提交内存耗尽而失败的。启动前先确认 `FreeVirtualMemory` 足够。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_stage.py run_stage0
& $py scripts\run_stage.py run_reps 1 1 50     # 场景 1，实现 1–50；按 worker 切分
& $py scripts\run_stage.py summarize
```
