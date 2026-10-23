# 复杂模型对比：从这里开始

一个2.94km宽、2.2km深、99道的复杂地质模型，两种噪声条件。模型含密集层系、褶皱、双断层、河道、薄层尖灭、低Q透镜体；观测含混合相位、非平稳衰减、多次波、相关噪声、相干干扰和坏道。
这是二维属性驱动的局部一维射线/卷积近似，不是实测或完整二维弹性波模拟。粗几何是构造并假定已知的先验。

| 条件 | 原分步法误差 | 优化方法误差 | 优化改善的区块 |
|---|---:|---:|---:|
| Moderate | 0.61893 | 0.56815 | 8/9 |
| Strong | 0.7886 | 0.70137 | 9/9 |

- index.html：离线完整图集。
- RESULTS_ZH.md：完整报告与局限。
- figures/complex_geology.png：原始地质属性。
- figures/Moderate_comparison.png、Strong_comparison.png：道集及误差。
- analysis/metrics.csv、block_metrics.csv：全部指标及逐区块结果。
- truth/model_truth.npz：地质属性、参考、多次波、噪声、静校正/增益等。
- core_p1/inputs/complex_input.json：求解器实际收到的观测与标定/公开几何，无真实Q。
- runs/：实际执行日志、冻结源码快照、MAT和JSON结果。

复现：已有完整观测与求解器，可直接运行 scripts/run_local.py --entry run_complex_model --label replay --timeout 5400 --run-pointer <新的run指针文件>。MATLAB路径在run_local.py中按本机配置。
重新生成模型用scripts/build_model.py（固定种子），绘图和评分用scripts/plot_and_report.py；Python需要numpy、matplotlib、Pillow。
原有算法文件保持冻结，没有根据本轮结果重新调参。数值收敛不等于恢复准确，深部/低能量区域和失败区块均保留。
