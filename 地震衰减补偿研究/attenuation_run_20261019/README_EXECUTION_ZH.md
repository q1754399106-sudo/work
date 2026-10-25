# 衰减补偿研究执行与交付

本目录用于执行用户交接包中的研究任务。原始压缩包、解压来源和历史证据保留，候选实现位于 `candidate/`。所有 MATLAB 工作连接用户已共享的 `ResearchMATLAB`；不启动或关闭实例。实际完成情况以各运行目录的 `RUN_METADATA.json`、连接回执和原始 CSV 为准。

交付入口是新建交付目录中的 `EFFECTS_ZH.md`（一页结果）、`effects.html`（完整图表与证据）、`RUN_METADATA.json`（实际范围及限制）及 `SOURCE_MANIFEST.json` / `FILE_MANIFEST.json`（哈希）。开发集、历史回放和人工评价夹具分别保留，不充当新测试效果。

## 执行顺序

1. P0 预检与历史重放；实际回执保存在 `preflight/`。
2. 候选数学检查、MATLAB 自检、开发运行。任何源码或配置的科学行为变化须先重新运行开发，再冻结。
3. 以 `candidate/experiment/freeze_experiment.py` 封存已完成的开发配置和源码，生成 `frozen_experiment_spec.json`。测试必须显式使用该文件，不读取测试结果调参。
4. 实际测试后，用 `candidate/audit/strict_evaluate.py` 核查完整矩阵、配对、模型划分、来源哈希、求解状态与预定工程门槛。负面结果同样完整报告。
5. 用 `candidate/visualization/render_results.py` 从实际测试 CSV、预定案例及保存的最坏案例绘图。图件的输入、变换、缺失值政策、描述文字和哈希写入 `FIGURE_MANIFEST.json`；视觉检查状态按实际清单报告。
6. 完整测试结束后运行 `audit_solver_logs.py`，生成 `experiments/<测试标签>_numerical_audit/`。它审查全部搜索求解的标志、正规方程残差、重启/回退、选中Direct/Split目标偏差以及参考/校准bootstrap边界，只生成后验文件，不改变冻结的CSV或元数据。原始CSV的 `TargetMismatch` 两路线均记录Direct目标偏差，Split偏差单独从其日志恢复。
7. 最后运行下列交付组装命令。参数均为目录；`--test-run` 指向包含 `metrics.csv` / `RUN_METADATA.json` 的 MATLAB 输出子目录，`--output` 必须不存在。

```powershell
$researchPython = 'C:\Users\qq175\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $researchPython .\build_delivery.py `
  --test-run .\experiments\test_v1\matlab `
  --evaluation .\experiments\test_v1_strict `
  --figures .\experiments\test_v1_figures `
  --numerical-audit .\experiments\test_v1_numerical_audit `
  --output .\delivery_v1
```

示例目录名仅说明接口，应使用实际运行路径。组装脚本只读取证据，不调用 MATLAB、不运行实验、不改变参数。它把输入、图件、执行源码、预检回执及冻结记录复制到新交付目录，生成在目录内可相对跳转的 HTML。

`--numerical-audit` 可以省略；此时自动查找测试运行同级的 `<测试标签>_numerical_audit`。交付器校验 `audit_checks.json` 中原始metrics及运行元数据的哈希，绑定成功后才展示诊断摘要。缺失或不匹配会明确标为未绑定，保留原审计文件但不据其作数值解释。后验目录中的诊断CSV、JSON、说明和审计工具源码均纳入交付；参考/bootstrap记录跨套件重复，不能解释为同等数量的独立校准或地质样本。

## 结果约束

- 三条候选：`NoiseUncertaintyOnly`、`QModelErrorOnly`、`CombinedUncertainty`；主路线为 Split，同时报告共享先验的 Direct。三组机制对照为 `ScalarBudgetDamping`、`NoiseCovarianceOnly`、`NoStructureCombined`。
- 候选 Split 求解完整非稳态映射 **M·Gref**。这相对历史分步实现有定义变化；不能假定非稳态算子的精确乘积分解。历史 Frozen **6.48%** 为原归档口径，不是本次测试指标。
- 候选 Direct 按 **Cn⁻¹** 对道加权，正则项采用 **λj/σref,j²**；继承的 Frozen Direct 使用原始数据范数及 **μ=γ·median(λ)**。因此相对 Frozen 的差异包含加权方式和路线定义变化，不能将其全部归因于 Q 误差协方差。`NoiseCovarianceOnly` 使用点估计 MAD/AR(1)，并非严格的白噪声点尺度单因素对照；应结合全部分支、Direct/Split 和机制对照解释。
- Main18 与 LowerSNR12 各用 `SmoothAVO`、`FaultSection`、`PinchoutSection`、`WeakSection` 四个白噪声场景。`CoherentNoiseAVO` 独立位于 `StressCoherent`，避免主套件无相干污染时重复计数。
- 四个测试模型是预选合成参数留出，尚无独立实测地质证据。Main/LowerSNR 每模型三个校准、每校准六个观测；压力套件仅一个校准、六个观测，重复层级较弱。
- bootstrap 百分位范围仅是当前模型/校准/观测设计下的条件性描述，不是总体置信区间或显著性检验。
- 任何暂定数值通过都不自动替换 Frozen。逐模型退化、局部保幅、最坏案例、失败行、差异原则端点与偏离目标记录均须查看。

## 组装校验与可移交性

`build_delivery.py` 将实际测试状态、严格评价的三项输入哈希、完整矩阵行数/观测数以及图件输入绑定后才允许产生改善表。实际完成但工程门槛失败是有效结果；缺失、未完成、人工夹具或哈希失配会生成明确的 `PARTIAL_NO_EFFECTIVENESS_CONCLUSION`，并写出原因，效果数值不补造。检查范围不等同于外部科学审计。

`SOURCE_MANIFEST.json` 逐文件记录实际执行的来源路径、原记录哈希和交付复制哈希；`FILE_MANIFEST.json` 与 `SHA256SUMS.txt` 记录交付文件。为避免递归哈希，它们不对自身求哈希；源码中的绝对历史路径属于溯源记录，重新执行应按当前目录重新冻结。压缩交付包时保留这些清单，排除 `dependencies/` 与缓存。

相关方法及假设见 `candidate/matlab/P1_THEORY_AND_ASSUMPTIONS.md`、`candidate/experiment/METHODS_ZH_EN.md`；文献差异与已核对公开来源见 `LITERATURE_DIFFERENCES_ZH.md`。报告附有状态为 `UNVERIFIED` 的 Material Passport，明确本地字节/数值核验与外部学术完整性核验的边界。
