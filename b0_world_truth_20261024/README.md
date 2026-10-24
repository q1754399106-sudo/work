# B0 的世界真值修正与未见模型确认

**状态：已暂停。** 开发阶段已完成，新模型 U1–U4 已生成；确认阶段没有运行，U1–U4 保留为未见数据。详见 [PAUSED_ZH.md](PAUSED_ZH.md)。`RESULTS_ZH.md`、`predictions/` 与 `analysis/truth_U*.csv` 目前都不存在。

开发声明：[DEV_PLAN_ZH.md](DEV_PLAN_ZH.md)；确认方案：[PLAN_ZH.md](PLAN_ZH.md)；哈希：`logs/DEV_PLAN_SHA256.txt`、`logs/PRE_CONFIRM_HASHES.json`、`predictions/SEALED_SHA256.json`。

## 内容

- `scripts/`
  - `diagnose2.m`：同一组噪声抽样下的三种世界真值：B0（FullCalRefit 输出）、B0M（自助均值）、B0C（第三轮 Standard 相干收缩）。
  - `coherence_reference.m` 与 `round3_20261024/run_combo.m` 的 refit_state、coherence_scores、apply_gain 逐句一致。
  - `make_unseen.m` 生成 U1–U4，正演与辅助函数逐字取自 `bae_q_20261024/make_independent.m`。
  - 其余辅助函数取自 `qscan_bias_diagnostic_20261024/scripts`。
- `inputs/`、`truth/`：新模型 U1–U4 的求解侧输入（公开配置 + 48 个诊断道集）与真值，分开保存。
- `analysis/`：开发结果（`dev_w*.csv`、`dev_selection.json`）、确认真值（`truth_U*.csv`）与评分。
- `predictions/`：确认预测与封存收据。
- `figures/`、`logs/`。
- `cache/`（不上传）：从 Git 对象逐字节取出的输入快照，以及每个道集的两个输出。

## 复现

需要 MATLAB R2024b、Python（NumPy、pandas、Matplotlib），以及仓库中相关目录的 Git/LFS 对象（`snapshot_inputs.py` 直接读对象，不依赖稀疏检出）。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_stage.py run_dev 1 6          # 2 6 … 6 6 并行
& $py scripts\analyze_dev.py                     # 按 DEV_PLAN 规则选定变体与阈值
& $py scripts\hash_preconfirm.py
& $py scripts\run_stage.py make_unseen          # 拒绝覆盖已有新模型
& $py scripts\run_stage.py run_confirm 1 6      # 2 6 … 6 6 并行
& $py scripts\seal_confirm.py
& $py scripts\run_stage.py run_truth 1          # 2、3、4
$env:PYTHONPATH='D:\ChatGPT\研究\attenuation_run_20261019\dependencies'; & $py scripts\analyze_confirm.py
```
