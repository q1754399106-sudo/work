# P1 独立合成实验 / Fresh synthetic experiment

这些文件是可运行的实验入口与冻结协议。本目录本身没有 MATLAB 反演结果；只有实际共享 ResearchMATLAB 会话写出的 `RUN_METADATA.json` 与全量 CSV 才可证明运行。

`p1_run_experiment(stage,outdir,specfile)` 只接受 `development` / `test`，只在新目录写输出。它保存并恢复 path 与 RNG，不使用 base workspace、不启动或停止 MATLAB，不创建或关闭MATLAB图窗。图在导出数据后由Python绘制。

使用现有共享 MATLAB 会话，先把本目录加入会话路径，然后调用：

```matlab
p1_run_experiment('development', 'C:/path/new_development_output', ...
 'C:/path/candidate/experiment/experiment_spec.json');
```

开发完成后，从 Python 生成冻结包。该命令核验开发运行、配置和源码一致，失败行不阻止冻结但会留在记录中；冻结不代表候选通过验收。

```text
python freeze_experiment.py --spec experiment_spec.json --development-run C:/path/new_development_output --output C:/path/new_freeze
```

冻结包内 `frozen_experiment_spec.json` 带显式 `freeze.expected_hash`、原始配置快照哈希、所有核心/候选源码哈希、开发记录哈希。只有这份文件能作为测试参数：

```matlab
p1_run_experiment('test','C:/path/new_test_output', ...
 'C:/path/new_freeze/frozen_experiment_spec.json');
```

测试前逐项核验配置与源码；源码、配置、种子、分组或阈值变化会使运行拒绝。发现测试结果以后不能以同一冻结设计调参。若确需改变方法，原测试应如实归为已揭盲，另选新的独立测试设计，保存原失败。

测试包含 4 个此前未使用的确定性合成分层参数集，主试验每模型 3 个独立校准种子，每校准 6 个观测种子，18/12 dB 各 4 个白噪声场景。8 个独立因素压力组各用每模型 1 校准×6观测×1预定场景；该部分不声称验证了校准方差。开发使用 2 个不同参数模型、18/12 dB、1校准×2观测×3场景，超参数搜索预算为零。所有方法都对同一观测输出 Direct 与 Split，保留普通阻尼与结构对照。固定矩阵为开发24观测/384行，测试768观测/12288行。

主/低信噪比组均只加白噪声。与原任务书每组列5场景的机器清单不同，冻结设计把 `CoherentNoiseAVO` 单列在 `StressCoherent`：因为它不加相干项时与相同种子的SmoothAVO完全重复，会重复加权同一观测。五种场景在完整协议中全部保留。所有污染均与名义白噪声 SNR 分开记载，不能把有色、相干、脉冲污染混入主试验均值。

`metrics_partial.csv` 每个方法后落盘；最后生成 `metrics.csv`。失败、NaN、迭代标志、正规方程残差、端点与目标失配均保留。完整求解器记录在 `solver_logs/`；估计Q/频率与边界在 `calibration/`。事先固定第一模型、第一校准、第一观测的示例、固定道1/6/11以MAT和JSON保存在 `examples/`。图由根任务在MATLAB外绘制，图题、坐标与图例应中英双语。`p1_plot_example.m`只是未调用的参考函数，主入口不使用它。

`RUN_METADATA.json` 含所有源文件哈希、配置文件原始字节哈希、层级种子、实际 MATLAB 版本和真值审计边界。图与CSV中的 `HTruthRatio`、`EstimatedSigmaRatio` 都是后验诊断，绝不回馈选参。

`worst_cases/` 按事先固定规则保存每个方法/路线最大有限RelativeError案例（包含数值仍有限的求解器失败），最多16个JSON/MAT案例；同时保存同一观测的Frozen Direct/Split用于比较。非有限失败的键和值保留在CSV与日志。这个后验挑选只服务最坏案例报告，不能反馈选参。

测试完整后用同目录的 `strict_evaluate.py` 做工程门槛和分层描述区间；不能以平均改善替代局部指标、失败率与 Direct 对照。外部地质证据仍缺失，任何门槛通过最多是 `PROVISIONAL_PASS`。
