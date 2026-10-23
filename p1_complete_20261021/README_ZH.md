# P1 完整重跑

本目录从用户此次提供的 ZIP 内 MATLAB 实现建立。旧 ZIP、旧结果及其他修订目录均保留。

计划：9 个模型、3 个校准种子、每次校准 6 个观测种子，七个方法、Direct/Split 两条路线。Main 与 12dB LowerSNR 各五个场景；30° 相位、+3% 速度、4ms 残余 NMO、有色噪声四种压力各两个 AVO 场景。合计 2,916 组观测、40,824 行指标。

速度压力条件为公开速度乘 1.03、层厚固定，重新校准 Q/主频并重建算子；不同于旧报告速度和层厚同时扰动的条件。残余 NMO 的填充 IFFT 结果先截取为观测长度，避免维度错误。

```powershell
py -3 scripts/run_complete.py --stage preflight --workers 1
py -3 scripts/run_complete.py --stage full --workers 2
py -3 analysis/analyze_run.py runs/<full_run>
```

每次运行输出到新的 runs 子目录。保存执行前计划、源码快照及哈希，每组观测完成后追加指标并更新进度。子进程独立启动，保留用户正在使用的 MATLAB。实验完成与验收通过分别记录。

实际覆盖、退出状态、未通过项及置信区间以生成的 EXECUTION_SUMMARY.json 和 analysis/RESULTS_ZH.md 为准。不得把代码存在或进程已启动写成实验已完成。

既有九个模型已被观察过，本次是完整交叉回归验证。即使运行计划和源码在本次执行前冻结，也不能将这些模型重新称为未接触盲测或独立外部证据。没有根据恢复真值在线挑选分支或剔除失败样本。
