# 失真预测诊断

先读 [RESULTS_ZH.md](RESULTS_ZH.md)。结论：仅反射系数敏感性有附加信息；完整 Q/主频/相位/校准指标未通过，不进入新补偿算法开发。

本目录保存 2026-10-24 已实际运行的一轮结果。`predictions/SEALED_SHA256.json` 绑定预测程序、白名单输入、冻结方案和预测文件；`analysis/` 为封存后真值评价。零观测窗的风险未定义，原始数值占位不得解释为可靠，详见 `analysis/prediction_domain.csv`。

可直接读取 CSV 和 PNG/PDF，无需 MATLAB。Python 分析需要 NumPy、pandas、Matplotlib；精确版本记录在报告中。已有环境位于仓库 `attenuation_run_20261019/dependencies`（按仓库规则不上传）。

```powershell
$env:PYTHONPATH='D:\ChatGPT\研究\attenuation_run_20261019\dependencies'
& 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'theory_predict_diagnostic_20261024\scripts\analyze.py'
& 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'theory_predict_diagnostic_20261024\scripts\write_report.py'
```

分析脚本校验封存文件，逐场景复现旧评分，生成对照表和六组图。需要恢复的旧表：`bae_q_20261024/analysis/dev_metrics_all.csv` 与 `ind_metrics_all.csv`。

完整计算需要 MATLAB R2024b，并恢复 `bae_q_20261024/source/`、`core/` 以及本目录 `inputs/` 的 Git LFS 文件。`scripts/run_stage.py` 可执行 `export_inputs`、`predict`、`score_outputs` 三个阶段，`scripts/check_calibration.m` 补充校准导数步长收敛检查。

请先把本目录复制到仓库内一个新的同级运行目录，再重算；移走新目录中的旧缓存、预测和分析产物，保留本目录原始封存。`export_inputs` 仅在需要重新从历史归档导出输入时运行，还需要恢复对应 `soft_refit_20261023`、`complex_fieldlike_20261023` 和 `bae_q_20261024` 原文件。它只复制白名单字段，不生成新模型。

本轮完整预测的读取真值前封存时刻记录在 JSON 中。历史标签已知，后续留一模型分析不能替代未见模型的最终验证。补充消融和绘图的后处理没有修改封存风险。
