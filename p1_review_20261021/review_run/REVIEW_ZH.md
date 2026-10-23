# P1 衰减补偿方法检查与试运行

检查日期：2026-10-21（Asia/Taipei）。原始 ZIP 与解压后的方法源码、历史结果均未修改；新增检查脚本和运行结果放在 `review_run/`。

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run
- Origin Date: 2026-10-21
- Verification Status: PARTIAL — 历史指标已重算；本机运行范围见下文，不代表全量复现
- Version Label: p1_review_20261021_v1

## 方法判断

方法有继续研究的价值，但目前应描述为“有效 Q 校准 + 数据驱动结构先验 + 经验模态正则化 + 相干干扰抑制”的合成数据实验。源码尚不足以支持“校准不确定性已经传播至模型误差协方差”的解释；交付脚本也存在妨碍完整复现的缺陷。

其基本流程为：从独立校准道集估计单一有效 Q 与主频；用观测道集的噪声窗 MAD 估计每道噪声尺度；通过保守初始反演学习结构算子；对 Offset 道集使用 AVO 正交投影，对 Section 数据使用局部尾端保护；以偏差原则选择共同尺度，分别求 Direct 与 Split 解。使用估计参数建立算子、保持两路线相同先验，以及保留局部振幅评价，均有合理的工程动机。

## 已确认的问题

| 优先级 | 发现与影响 | 源码位置 |
|---|---|---|
| 高 | MATLAB 全量入口把 `uq_case` 的第三返回值当成 scene；第三值实际是 map，scene 是第四值。随后 `uq_observation` 对数值数组访问 `.truth` 会报错。检查用驱动最初复用了这段调用，出现相同错误；后续仅修正检查驱动。 | `core_p1/run_p1_research.m:27`；`core_p1/matlab_src/uq_case.m:1` |
| 高 | 多次校准更新 `cm_c.Q/fpeak` 后，仍从 `invmodels_base` 取旧算子，未按新 Q 重建传播模型。因此重复校准结果没有完整传播到反演；不能据此充分证明校准稳健性。历史扩展和新地质循环均有此问题。 | `scripts/run_p1_scipy_evaluation.py:213–223, 294–305` |
| 高 | `StressVelocityError` 未修改 vp，也未重建传播算子；实际仅执行 `cm_st['fpeak'] *= 0.96`。不能把该结果称为速度误差测试。配置声明 +3% 速度、中文报告声称 4%，执行代码均不相符。 | `scripts/run_p1_scipy_evaluation.py:230–238` |
| 高 | `p1_calibrate_ensemble` 计算 `log_q_std/log_f_std/log_qf_corr`，但 `uq_inverse_config` 只保留 Q 和主频点估计；求解器没有使用这些不确定性。Q 权重依赖奇异值、模态频率和手设常数；`b_eD/b_eS` 始终为零。函数注释所称 block-bootstrap ensemble 也未实际执行。 | `core_p1/matlab_src/p1_calibrate_ensemble.m`、`uq_inverse_config.m`、`p1_q_model_error.m` |
| 高 | 无结构 MATLAB 分支的初始 SVD 解没有包含 `extra_diag`，但正规方程残差却包含它。Section 尾端保护非零时，返回解不再求解所声明的方程。Python 的无结构对照直接传入 `None`，两版本在此也不一致。 | `core_p1/matlab_src/p1_solve_branch.m:96–112`；`scripts/p1_solver_core.py:388–394` |
| 中 | `QModelErrorOnly` 仍会计算并启用 Section 的尾端噪声保护，而非单独开启 Q 误差项。分支消融不完全纯净。 | `core_p1/matlab_src/p1_solve_branch.m:28, 43–45, 98` |
| 中 | Python runner 与核心依赖未随包交付的 `verify_numpy_solver`，并硬编码 `/tmp/work/...` 和 `/working_dir/...`；原脚本无法从该 ZIP 独立启动。报告中的 `--action smoke/full` 也不在实际 CLI 的合法动作中。 | `scripts/run_p1_scipy_evaluation.py:8–11`；`scripts/p1_solver_core.py:12–13`；`scripts/run_matlab_p1.py:9` |
| 中 | 新地质每个校准只有每场景一个观测种子；历史额外两次校准也只有每场景一个种子，未达到声明的每校准至少 6 次观测。新地质的校准种子与观测种子同步变化，难以分离两者影响。 | `scripts/run_p1_scipy_evaluation.py:220, 301`；本次 `independent_audit.json` |
| 中 | `StrongerDamping` 只比较提高目标到 1.15，而候选常使用低于 1 的目标，同时改变频谱权重和结构强度。此消融不足以“严格证伪”更一般的阻尼/频带调整解释，应增加相同目标与相同先验的对照。 | `core_p1/inputs/p1_experiment_spec.json` |

真值隔离方面，反演入口确有 public 字段白名单，未在已读的反演函数中看到直接读取 `scene.truth`。但包内 Python “动态投毒审计”只比较本地重新写出的字段筛选函数，没有重新执行真实求解器，也未扫描 Python 核心；不能将其提升为完整动态端到端验证。公开输入仍假定层速度和厚度可用，这属于实验条件，不能外推为完全盲反演。

## 历史结果重算

通过包内 `evaluate_candidate.py` 对原 CSV 重新计算，未重新生成这些历史波形：

| 范围 | 联合分支结果 |
|---|---|
| 历史集 126 对 | `PROVISIONAL_PASS`，39/39 项工程检查；主试验等模型误差从 0.685793 降至 0.635630，改善 7.3147% |
| 扩展集 246 对 | `FAIL`，113/115 项通过；误差 > 1 的次数从 Frozen 的 17/246 降为 0/246 |
| 局部未通过 1 | `MultiCalMain/SoftLowQ/AVOGradientError`：0.620129 → 0.654486，增加约 5.54%，超过 5% 容许值 |
| 局部未通过 2 | `NewGeologyMain/Holdout_ThinInterbed/WeakEventError`：1.351831 → 1.439070，增加约 6.45%，超过 5% 容许值 |

扩展集报告中的 7.3147% 仍是 `Suite == Main` 的历史三模型统计，不是九模型全集平均改善。工程门槛通过不等于统计显著，零失效仅表示这批有限样本没有超过指定阈值。

交付清单列出的 136 个文件均通过大小与 SHA-256 核对。原说明提及的另一批 156 文件、36 个 MATLAB 历史例子的独立数值双生验证依赖未交付内容，本次没有将历史声明当成重新验证的事实。

## 本机试运行

实际环境：Windows，MATLAB R2024b Update 6；辅助检查使用 Python 3.12.14。现有 MATLAB 未发布共享 `ResearchMATLAB`，因此包内共享会话入口不能连接。使用独立隐藏后台 MATLAB 批处理执行检查，1200 秒硬超时，输出另存。

包内 `validate_p1_unknown_q` 已实际通过：公开输入字段投毒一致性、同模板 Q/主频恢复，以及 Frozen/Combined 的正规方程残差。该投毒测试验证公开输入构造相同，并未在两个不同私有模型上重复整个补偿求解；同模板恢复也属于模型内一致性测试。

检查驱动覆盖：SoftLowQ 一个模型、校准种子 20261030、观测种子 20261101、18 dB、五个场景、七个方法、Direct/Split 两路线。完整试跑状态、逐场景结果及数值残差附于下方。该小样本的联合分支通过包内工程验收，但不能代替全量评估。

运行命令：

```powershell
& 'C:/Users/qq175/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' review_run/launch_matlab.py
```

首次驱动中调用 `uq_case` 的返回值错误被记录在 `matlab_stderr.log`；修正后的驱动使用第四返回值，输出到 `matlab_smoke_corrected/`。原包源码没有被修复。本次没有执行九模型全量网格，也没有宣称 Python 与 MATLAB 的 P1 全分支严格一致；相同数字种子在两种随机数实现下也不保证相同噪声样本。

## 建议处理顺序

1. 修正 MATLAB 返回值和无结构分支；补齐 Python 核心依赖并改用包内相对路径。
2. 每次校准后重建算子和 AVO 投影；速度压力测试实际扰动速度并重建算子。
3. 明确选择“经验模态正则化”定位，或真正传播 Q/主频不确定性；完善纯 Q、纯噪声、相同目标的消融。
4. 冻结上述实现后，按校准 × 观测种子完整交叉重跑独立模型，重新计算局部保幅和分层区间。

## 本机试跑实测结果

状态：COMPLETED；70 行指标；MATLAB 内计算耗时 178.1 秒。

| 场景 | Frozen Split 相对误差 | Combined Split 相对误差 | 改善 |
|---|---:|---:|---:|
| CoherentNoiseAVO | 1.132113 | 0.779819 | +31.12% |
| FaultSection | 0.666569 | 0.669550 | -0.45% |
| PinchoutSection | 0.625761 | 0.611701 | +2.25% |
| SmoothAVO | 0.667976 | 0.644100 | +3.57% |
| WeakSection | 0.660892 | 0.666129 | -0.79% |

五场景平均误差 0.750662 → 0.674260，相对改善 10.18%。仅是单模型单种子的描述性结果。

联合分支两路线最大正规方程相对残差 9.564e-08。超过 2e-7 门槛或求解标志非零的记录共 6 行：

- FaultSection / CombinedUncertainty_NoStructure / Direct：NormalResidual=1.990e-04，FinalFlag=0。
- FaultSection / CombinedUncertainty_NoStructure / Split：NormalResidual=1.970e-04，FinalFlag=0。
- PinchoutSection / CombinedUncertainty_NoStructure / Direct：NormalResidual=2.259e-04，FinalFlag=0。
- PinchoutSection / CombinedUncertainty_NoStructure / Split：NormalResidual=2.244e-04，FinalFlag=0。
- WeakSection / CombinedUncertainty_NoStructure / Direct：NormalResidual=7.099e-05，FinalFlag=0。
- WeakSection / CombinedUncertainty_NoStructure / Split：NormalResidual=6.950e-05，FinalFlag=0。

无结构分支的实际残差失败与源码缺陷一致；不能因为进程 exit code=0 或 FinalFlag=0 就判全部求解通过。详见 `smoke_summary.json`、`matlab_smoke_corrected/metrics.csv` 与 `matlab_run2_stdout.log`。
