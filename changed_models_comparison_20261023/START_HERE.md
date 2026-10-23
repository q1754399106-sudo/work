# 本轮更换模型对比

已完成4个新地质配置、80组普通噪声及20组有色噪声，原有求解代码保持冻结。

普通噪声整体误差：原分步法0.8377，优化方法0.1755；优化方法在76/80个案例更低。
普通噪声可辨识弱反射平均误差：1.6704 → 0.5011；1个模型弱反射退化。
有色噪声弱反射平均误差：1.8184 → 0.7653。

- RESULTS_ZH.md：完整报告，包括失败案例和适用范围。
- index.html：离线道集图库，全部100组均可查看。
- figures/comparison.png：逐模型误差对比。
- figures/weak_amplitudes.png：4个固定案例的弱反射振幅曲线。
- figures/worst_relative_case.png：按误差比自动选出的最差案例。
- metrics.csv：原始500行评分；metrics_identifiable.csv：统一排除不可辨识弱振幅后的评分；summary.json：两种评分及汇总。
- PLAN.json / FROZEN_SOLVER_HASHES.json：运行前登记计划与原有求解代码哈希。
- runs/：MATLAB日志、源码快照、独立校准和原始MAT/JSON。

主对照是原分步法FixedSplit；FixedSplitSharedCal另使用与优化方法相同的完整校准结果。
四层、11道、Ricker子波正演体系保持不变；本轮检验新地质配置，不代表实测资料验证。
薄层WeakSection中央道有两个反射完全重合。弱反射评分在运行后复核时修正，统一排除这一不可辨识分量；原始评分保留，完整波形评分不变。详见报告的评分修正说明。

新增：avo_comparison/index.html 为40组保留/去除AVO对比图库，处理定义见 avo_comparison/METHOD_ZH.md。
