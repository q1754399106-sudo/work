# QScan 偏差—方差来源分解

先读 [RESULTS_ZH.md](RESULTS_ZH.md)。冻结方案：[PLAN_ZH.md](PLAN_ZH.md)；冻结后改动：[logs/POST_FREEZE_CHANGES.md](logs/POST_FREEZE_CHANGES.md)。

## 内容

- `scripts/`：MATLAB 阶段脚本与 Python 分析。求解器、噪声生成和 σ 估计全部调用 `bae_q_20261024` 原函数；`qscan_select.m` 是 `qscan_refit.m` 选择步骤的逐句拷贝，另外保留目标函数两项，阶段 0 验证其与原函数逐位一致。
- `analysis/`：全部数值表。`decomposition.csv`（每个估计量与每个固定因子的 MSE/偏差²/方差，因子扫描也在其中）、`qscan_vs_base.csv`（配对差）、`window_moments.csv`（35 道 × 26 窗）、`window_joined.csv`（再并入封存的 U 指标）、`time_profile.csv`、`spectral_transfer.csv`、`residual_attenuation.csv`、`replicates.csv`（每个实现的误差与所选因子）、`snr_probe.csv` 与 `snr_trend.csv`（补充登记 A）、`decisions.json`（D1–D5 与 P-A）。
- `results/`：`rep_*.mat` 为每个场景 100 个实现的逐样本一阶/二阶累加量和逐实现标量；`clean_*.mat` 为无噪声对照输出；`snr_*.mat` 为高信噪比探针。逐样本均值/方差可由 `summarize.m` 从 `rep_*.mat` 精确重建，存于本地 `cache/`，不上传。
- `figures/`：8 组图，各有 PNG、PDF 与说明文本。
- `logs/`：各阶段原始输出、执行记录、冻结哈希。

## 复现

需要 MATLAB R2024b，以及仓库中的 `bae_q_20261024/{core,source,inputs/independent,truth/independent,results/ind,results/ind_qscan}`（稀疏检出需先恢复）。Python 分析需要 NumPy、pandas、Matplotlib（仓库 `attenuation_run_20261019/dependencies`，按规则不上传）。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_upstream.py        # 可选：把上游文件按哈希拷到 cache/，避免运行中稀疏检出变化
& $py scripts\run_stage.py check_reproduction
& $py scripts\run_stage.py run_clean
& $py scripts\run_stage.py run_replicates 1 1 100   # 场景 1–4 可并行
& $py scripts\run_stage.py summarize
& $py scripts\run_stage.py run_snr_probe 2 20 1 40   # 补充登记 A：模型号 1/2 × 20/26 dB
& $py scripts\run_stage.py summarize_snr
$env:PYTHONPATH='D:\ChatGPT\研究\attenuation_run_20261019\dependencies'; & $py scripts\analyze.py
```

`run_replicates` 拒绝覆盖已有分块结果；重算请在新目录进行，保留本轮记录。
