# 联合正则重构实验结果

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run
- Origin Date: 2026-10-22
- Verification Status: ANALYZED; numerical invariants checked; no independent replication
- Version Label: joint_refactor_v1

完成 640 行，MATLAB 退出码 0，当前数值批次耗时 1094.4 秒。四开发模型、同一既有校准种子、两个新观测种子；所有结果均为开发诊断。

本轮新校准分步候选均未通过预定工程筛选，不提升生产默认。

本批数据中，保留原绑定规则时，新校准使平均误差降低 10.11%；在新校准上进一步归一化并解耦，较小固定结构系数版只再改善 0.75%。收益主要体现为校准替换的全流程效应，尚不足以支持推广新正则。

## 新校准下的分步结果

| 方法 | 平均相对误差 | 相对 Legacy 改善 | >1计数 | 最差 | 端点数 |
|---|---:|---:|---:|---:|---:|
|Legacy|0.610439|+0.00%|0/32|0.831679|0|
|DecoupledWeak|0.605887|+0.75%|0/32|0.831869|0|
|Decoupled|0.614887|-0.73%|0/32|0.872188|3|
|SharedIsotropic|0.612525|-0.34%|0/32|0.835819|0|
|Physical|0.684303|-12.10%|0/32|0.923691|12|

## 校准与正则的交叉结果（分步）

| 方法 | 旧校准误差 | 新校准误差 | 校准净改善 |
|---|---:|---:|---:|
|Legacy|0.679095|0.610439|+10.11%|
|DecoupledWeak|0.675643|0.605887|+10.32%|
|Decoupled|0.682138|0.614887|+9.86%|
|SharedIsotropic|0.678161|0.612525|+9.68%|
|Physical|0.665417|0.684303|-2.84%|

## 同共享正则对照的物理惩罚净作用

| 校准 | 路线 | SharedIsotropic | Physical | 加入RQ的相对改善 |
|---|---|---:|---:|---:|
|Old|Direct|0.678161|0.665415|+1.88%|
|Old|Split|0.678161|0.665417|+1.88%|
|New|Direct|0.612525|0.684283|-11.72%|
|New|Split|0.612525|0.684303|-11.72%|

## 新校准下局部指标（分步，越低越好）

| 方法 | AVO截距 | AVO梯度 | 弱反射 | 尖灭尾部 |
|---|---:|---:|---:|---:|
|Legacy|0.373460|0.407020|0.893851|0.002859|
|DecoupledWeak|0.362969|0.396858|0.904526|0.003054|
|Decoupled|0.342621|0.375503|0.934357|0.007241|
|SharedIsotropic|0.381404|0.409281|0.890383|0.003052|
|Physical|0.316659|0.313872|0.881333|0.004909|

## 工程筛选与数值检查

- Direct / DecoupledWeak: FAIL，整体改善 +0.75%，失败条件 4 项。
- Direct / Decoupled: FAIL，整体改善 -0.71%，失败条件 8 项。
- Direct / SharedIsotropic: FAIL，整体改善 -0.31%，失败条件 10 项。
- Direct / Physical: FAIL，整体改善 -12.06%，失败条件 11 项。
- Split / DecoupledWeak: FAIL，整体改善 +0.75%，失败条件 4 项。
- Split / Decoupled: FAIL，整体改善 -0.73%，失败条件 8 项。
- Split / SharedIsotropic: FAIL，整体改善 -0.34%，失败条件 9 项。
- Split / Physical: FAIL，整体改善 -12.10%，失败条件 11 项。

全部搜索候选最终求解标志失败数 0；最大法方程相对残差 9.999e-08；旧绑定/归一化重放差 4.035e-05；期望失配恒等式误差 1.209e-15；共享正则矩阵最大相对差 0（逐列实际计算）；运行快照变化 0 个文件。

## 解释边界

归一化+解耦相对 Legacy 是组合消融，未分别在全网格运行仅归一化与不归一化解耦；仅归一化且补偿系数的等价性通过重放检查。相干预减、局部 E、先导和 H 保持既有实现。SharedIsotropic 与 Physical 使用相同 Gbar，因此它们的差异可以归于 RQ；Physical 与 Legacy 的差异同时包含模态形状、数据算子和参数解耦。
旧/新校准使用相同正则构造规则；校准变化后，先导、lambda、H、经验模态权重和残差预算随规则更新。因此校准交叉反映校准替换的全流程效应，不是固定数值正则矩阵只替换 G 的消融。
Legacy 是相同 CalShrink025Light 配置下的原绑定求解，不是最早 Combined 或 Frozen；本轮没有重新运行 Frozen 全网格比较。DecoupledWeak 的 0.3 仅表示其归一化结构系数小于 Decoupled 的 1，并不保证它比原绑定方案的实际结构惩罚更弱。
RQ 的情景范围来自未校准的不确定性代理。Split 固定 K、使用 K Gbar 数据项，但采用共同的未白化 RQ；这是控制实验定义，不是 Split 白化度量下的期望风险。未测试误差协方差似然、稳健结构、速度压力、新校准种子、全流程 AVO 保护和实测资料。
新观测种子不能把已用于开发的模型和校准变为独立确认。两种子描述性均值不支持显著性或后验覆盖率主张。alpha 由 Direct 残差选定后共享给 Split，不能据此判断各路线最优能力。未调整原确认门槛，未提升生产默认。

逐行指标：C:\Users\qq175\Documents\ChatGPT\研究\joint_refactor_20261022\runs\paired_continuation_20261022T113517Z\matlab\metrics.csv
设计：PLAN_ZH.md；聚合：analysis/summary.csv；逐模型：analysis/model_summary.csv；全部局部失败：analysis/local_gate_failures.csv；全部筛选原因：analysis/gates.json。

物理惩罚各项数值及残差失配见 analysis/physical_penalty_diagnostics.csv。惩罚项的数值比例用于诊断，不能直接解释为独立的误差贡献或因果效应。

独立稠密解核对：{"completed": true, "max_dense_oracle_relative_error": 3.0703238157220953e-07, "routes_tested": 2, "methods_tested": 4}；证据 C:\Users\qq175\Documents\ChatGPT\研究\joint_refactor_20261022\runs\algebra_check_20261022T103808Z\matlab\checks.json

研究入口端到端核对：{"completed": true, "finite_output": true, "truth_field_rejected": true, "solver_flags": 0, "max_normal_residual": 9.775108696790344e-08, "shared_regularizer_difference": 0}；证据 C:\Users\qq175\Documents\ChatGPT\研究\joint_refactor_20261022\runs\api_check_20261022T105102Z\matlab\checks.json

运行批次：[{"run": "paired_continuation_20261022T113517Z", "status": "COMPLETED", "exit_code": 0, "elapsed_seconds": 1116.0160000000033}, {"run": "paired_development_20261022T103456Z", "status": "HARD_TIMEOUT", "exit_code": null, "elapsed_seconds": 3600.452999999994}]。
首批触及预设一小时时限后，保留完整观测组检查点并补跑缺失组合；完成行直接保留，不是一次不中断的640行运行。续跑数值求解与物理算子源码逐字节核对相同。上述最大法方程残差包括保留行与新行。
保留行逐字段核对：[{"source": "paired_development_20261022T103456Z", "retained_rows": 500, "max_numeric_relative_difference": 0}]。

图件与逐点数据：figures/calibration_regularization.png、figures/local_paired_ratios.png、figures/model_plot_data.csv。图件各有 PDF 和说明，来源与图像哈希见 figures/manifest.json。
