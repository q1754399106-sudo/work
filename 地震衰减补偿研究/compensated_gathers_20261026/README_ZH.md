# 当前推荐方法的补偿道集

手机入口：[在线图册](README.md)。覆盖前一册的 29 个模型配置，另有六个 Q 低估 30% 的 H 模型场景。每个场景固定一个 14 dB 噪声实现，不能把这些展示当成新一轮多重复统计优势检验。

14 dB 为名义随机噪声设定；V3/V5 另有相干噪声、多次波、静校正扰动及脉冲干扰，总干扰信噪比不等于 14 dB。观测直接沿用归档，不删去这些干扰。

## 方法与含义

按项目已有推荐规则绘制：Q 可信时 SC，Q 不确定时 QScan＋SC。这里的 SC 是 novelty_evidence_20261025 中对文献目标的 FISTA 独立实现，不是论文作者的原程序。SC 固定 a=0.125、mu_relative=10、最多 5000 次迭代、KKT 容差 0.02；本次没有调参。F1/F2 是原 SC 选参数据，不能视为独立测试。

这套使用规则不等于“所有模型上都由同一方法取得最低真值误差”。V4/V5 原研究已经有自研方法更好的场景，本册保留 V3/V4/V5 自研 Standard / QScan＋Standard 的道集对照，并显示实际误差，不隐去这种差异。主图的方法在计算前由 PLAN_ZH.md 固定，真值仅用于评分与绘图。

## 数据

- 四层十配置：原 make_stability_scene('FaultSection',c)，11 道，横向间距 50 m，偏移距 750 m；原 ray_absorption / attenuation_matrix 正演。固定种子 2026102601 至 2026102610 加入目标噪声水平 14 dB 的白噪声；噪声标准差使用生成时已知值。这些是本次新生成的展示观测，模型参数与前一图册一致。
- V1/V2/V4：原 Moderate 观测、给定 Q 的 SC。V4 同时提供原 MATLAB Standard。
- V3/V5：原 Moderate 观测、四层标定 Q 的 QScan＋SC。V3 由原 q_compensate 的 QScan 规则新计算因子；V5 复用 round3 的 scan.factor 和 scanned.Q。没有使用细层真 Q 或目标挑选因子。
- F1–F4：原 Moderate 观测和已归档 SC 输出。
- U1–U4：原 SNR14、replicate=1 观测，新算 SC。
- H1–H6：sc_neutral 的原 L1、replicate=1 含噪道集；参考 Q 场景复用 SC，低估 Q 30% 场景复用 qscan_sc 的 QScan＋SC。两张图共享同一份含噪观测，改变的是补偿时使用的 Q 及 QScan 校正。

复用的输入/输出文件记录 SHA-256，新算结果记录完整求解信息；SOURCE_MANIFEST.json 保存固定源提交及源文件校验。results/*.mat 内有观测 D、无衰减目标 target、补偿 compensated、时间 t、噪声 sigma 和实际求解配置 config，可独立重画全部图，无需重新求解。

## 图与指标

四幅道集图共用对称色标，范围覆盖目标、观测、补偿与残差全部数值。没有 AGC、逐道归一化、分位数裁幅或后期增强。下方两条叠加波形均为原始幅度；实线为目标，虚线为补偿。

相对波形误差为 ‖X−T‖F / ‖T‖F；观测误差以 D 代替 X。保留全部时样和道。V3/V5 死道在观测和输出中为零，但目标中仍有响应，因此全道指标包含死道损失；另提供 output_live_error，不用活道指标冒充原先全道口径。

SC 迭代信息按原容差保留，未收敛情况如发生会在图中明确显示，不用改变阈值将其标为成功。

## 复现绘图

MATLAB R2024b：

```matlab
addpath('D:\ChatGPT\研究\地震衰减补偿研究\compensated_gathers_20261026');
plot_gathers
```

绘图程序对已有完整 PNG/PDF/FIG 组合跳过导出，缺失图会补画。若要改版重新导出，请先将需更新的旧图移出 figures/。PNG 为 180 dpi，PDF 使用矢量线条/文字与矩阵图像，FIG 可在 MATLAB 中继续编辑。

求解脚本 run_cases.py 依赖 NumPy 2.3.5、SciPy 1.18.1；读取与原研究相同的 BLAS 单线程求解器。快照脚本会从固定 Git 提交恢复上游输入。H 与 F 的复用分支需要本机既有求解缓存；即使没有这些缓存，已归档 results/ 仍足以重画全部交付图。
