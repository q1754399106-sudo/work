# QScan 选择偏移的模型间稳定性

**状态：已完成。** 先读 [RESULTS_ZH.md](RESULTS_ZH.md)：偏移在 −4 到 +4 个网格步长之间随模型变化，常数去偏不可行；无噪声最优在 8/8 个模型上都在上尺度 Q0。冻结方案见 [PLAN_ZH.md](PLAN_ZH.md)，哈希见 `logs/PRE_RUN_HASHES.json`。

上游：`combo_bias_variance_20261024` 阶段 2（F1/F2 上 QScan 比组合最优 Q 高约 2 个网格步长）。本轮在 8 个新开发模型 G1–G8 上检验这个偏移是否稳定，判断“常数去偏”是否可行。不使用 U1–U4。

## 内容

- `scripts/os_setup.m`：路径、常数、G1–G8 的种子与宏观模型、3 个场景、噪声种子规则。
- `scripts/make_dev.m`：生成 G1–G8。正演与辅助函数逐字取自 `bae_q_20261024/make_independent.m`（经 `b0_world_truth_20261024/scripts/make_unseen.m`）。公开配置与无噪声道集存于 `models/`，评分用的目标存于 `truth/`。
- `scripts/run_offset.m`：扫描 worker，与 `combo_bias_variance_20261024/scripts/run_sweep.m` 逐句相同（去掉阶段 1 核对）。每个实现单独保存，已存在的会跳过。
- `scripts/summarize_offset.m`：误差曲线、D/S1/S2/B/最优、偏移、无噪声最优、振幅增益，以及 O1–O3 判定。
- `scripts/make_figure.py`：`figures/offset_curves.png/pdf`。
- `scripts/qscan_select.m`、`noise_replica.m`、`estimate_sigmas.m`：与 `combo_bias_variance_20261024/scripts` 中的文件逐字相同。
- `scripts/snapshot_inputs.py`、`hash_prerun.py`、`run_stage.py`：输入快照、冻结哈希、带日志与超时的 MATLAB 运行器。
- `cache/`（不上传）：输入快照与每个实现的扫描输出（约 4 MB/实现）。

## 运行

需要 MATLAB R2024b 和 Python。每个 MATLAB `-nojvm` 进程约需 2 GB 内存，最多同时跑 2 个。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_stage.py make_dev            # 拒绝覆盖已有模型
& $py scripts\run_stage.py run_offset 1 1 0 6   # 模型 1、场景 1、实现 0–6；模型 1–8 × 场景 1–3
& $py scripts\run_stage.py summarize_offset
$env:PYTHONPATH='D:\ChatGPT\研究\attenuation_run_20261019\dependencies'; & $py scripts\make_figure.py
```
