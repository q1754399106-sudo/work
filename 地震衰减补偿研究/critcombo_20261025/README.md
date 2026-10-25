# 已有选 Q 准则在组合层面的重评

**状态：已完成。** 先读 [RESULTS_ZH.md](RESULTS_ZH.md)：没有一个准则通过 W1，已有准则在组合层面也不能替代 QScan。冻结方案见 [PLAN_ZH.md](PLAN_ZH.md)（提交 eb27e94，哈希 `logs/PRE_RUN_HASHES.json`）。

把 A0–A3（`qscan_criterion_20261024`）与 C1–C4（`qcriteria_20261024`）的选择接到相干标准版上，按最终误差与 QScan（A0）比较。数据：F1、F2、G1–G8 的 30 个场景 × 6 个实现，复用 `combo_bias_variance_20261024` 与 `qscan_offset_20261025` 的扫描缓存。不使用 U1–U4。

## 内容

- `scripts/cc_setup.m`：路径、常数、30 个场景与缓存位置；`load_model.m`：读取公开配置、无噪声道集，评分时读取目标。
- `scripts/run_crit.m`：worker。按缓存中的种子重新生成观测，在 13 个网格因子上计算各准则；A0 必须与缓存中的 QScan 选择相同。
- `scripts/qscan_stats.m`、`select_factor.m`：逐字取自 `qscan_criterion_20261024/core`。
- `scripts/criteria.m`、`spectral_factor.m`：逐字取自 `qcriteria_20261024/run_criteria.m`。
- `scripts/noise_replica.m`：与前几轮逐字相同。
- `scripts/score_crit.m`：读取缓存中的组合输出，计算 MSE、相对 A0 的变化、选择误差与 W1 判定。
- `scripts/make_figure.py`：`figures/criteria_combo.png/pdf`。
- `scripts/snapshot_inputs.py`、`hash_prerun.py`、`run_stage.py`：输入快照（含扫描缓存哈希）、冻结哈希、MATLAB 运行器。
- `cache/`（不上传）：输入快照与每个道集的准则结果。

## 运行

需要 MATLAB R2024b、Python，以及两个扫描缓存（`combo_bias_variance_20261024/cache/sweep`、`qscan_offset_20261025/cache/sweep`，可按各自 README 重建）。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py
& $py scripts\run_stage.py run_crit 1 30
& $py scripts\run_stage.py score_crit
$env:PYTHONPATH='D:\ChatGPT\研究\attenuation_run_20261019\dependencies'; & $py scripts\make_figure.py
```
