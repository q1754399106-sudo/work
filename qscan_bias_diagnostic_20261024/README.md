# 无真值偏差诊断

先读 [RESULTS_ZH.md](RESULTS_ZH.md)。冻结方案：[PLAN_ZH.md](PLAN_ZH.md)；运行前哈希：`logs/PRE_RUN_HASHES.json`；预测封存：`predictions/SEALED_SHA256.json`。

## 内容

- `scripts/`：MATLAB 阶段脚本与 Python 分析。
  - `diagnose.m` 是只读观测与公开配置的诊断核心（B0 自助、G 门限、CV 留出预测、RS 残差）。
  - `full_refit.m`、`qscan_select.m`、`noise_replica.m`、`estimate_sigmas.m` 逐字取自 `qscan_bias_variance_20261024/scripts`。
  - 求解器本体来自 `bae_q_20261024`。
- `predictions/`：`manifest.csv`（68 个道集）、各工作进程的诊断结果、封存收据。
- `analysis/`：输入核验、开发探针、真值重复实现、逐道集评分、准确率汇总、`decisions.json`。
- `figures/`：PNG、PDF 与说明文本。
- `logs/`：各阶段输出与执行记录。
- `cache/`（不上传）：`snapshot_inputs.py` 从 Git 对象逐字节取出的全部输入，以及每个道集的两个输出。它们都可由脚本重建。

## 复现

需要 MATLAB R2024b、Git LFS 中的 `qscan_criterion_20261024` 与 `round3_20261024` 数据对象，以及 Python（NumPy、pandas、Matplotlib）。

```powershell
$py='C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $py scripts\snapshot_inputs.py                 # 从 Git/LFS 对象建立 cache/inputs 并校验
& $py scripts\run_stage.py check_inputs
& $py scripts\make_manifest.py
& $py scripts\run_stage.py run_diagnostics 1 4    # 2 4、3 4、4 4 可并行
& $py scripts\seal_predictions.py
& $py scripts\run_stage.py run_truth_replicates qc F3 Moderate   # 8 个确认场景
& $py scripts\run_stage.py score_diagnostics
$env:PYTHONPATH='D:\ChatGPT\研究\attenuation_run_20261019\dependencies'; & $py scripts\analyze.py
```

这台机器被其他会话的 MATLAB 作业占满时，单线程（`bd_setup.m` 中 `P.threads=1`）比四线程快约三倍。
