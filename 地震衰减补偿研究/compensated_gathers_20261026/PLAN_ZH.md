# MATLAB 补偿道集图册：展示规则

任务：为上一册 29 个主要模型配置生成补偿道集，并发布手机可直接浏览的 GitHub Markdown 图册。

“最佳”按项目当前推荐的使用规则解释，不在看到真值后逐图选择最小误差方法，不把推荐方案称为所有模型的已证最优。

- OriginalControl 与 P1 九个四层配置：用归档 make_stability_scene 的 FaultSection（11 道，共同偏移距 750 m），准确 Q，SC；没有对应的现成观测，使用原正演算子生成观测，固定 14 dB 白噪声与种子 2026102600+配置序号。这些是新生成的展示场景，不冒充已归档实验。
- V1/V2/V4：原 Moderate 14 dB 观测，给定 Q 的 SC。V4 补充原 Standard 对照，已有证据显示其在该生成器上更好。
- V3/V5：原 Moderate 14 dB 观测，四层标定 Q 属于不确定输入，使用 QScan＋SC；V5 同时保留原 QScan＋Standard、未叠加 SC 等对照。V3 使用原 QScan 规则新计算公共 Q 因子；V5 复用 round3 已归档因子，不用真值挑 Q。
- F1–F4：原 Moderate 14 dB；U1–U4：原 SNR14、replicate=1；已知宏观参考 Q，使用 SC。F1/F2 是 SC 选参数据，图中和说明中标注。
- H1–H6：sc_neutral 原 Q100_14dB、replicate=1 的 SC 输出；另画原 Q070_14dB、replicate=1 的 QScan＋SC，展示 Q 低估 30% 的使用条件。
- SC 固定 a=0.125、mu_relative=10、KKT<=0.02、最多5000次，代码逐字取自 novelty_evidence 的 benchmark/common。优先复用已有输出，并校验文件哈希；缺失场景仅重算，不调参。求解与真值评估分离。
- 每张主图显示目标、观测、补偿、补偿减目标误差；共用对称幅度色标，标出完整幅度范围，不做 AGC、逐道归一化或分位数裁幅。另提供目标/补偿波形叠加，帮助看到弱反射和深层变化。
- 每个场景显示相对 Frobenius 波形误差，单个固定噪声实现，仅作展示，不新增统计优势结论。死道随原数据保留，图中明确指出。

导出 PNG、PDF、可编辑 FIG，生成 GitHub README 图片画廊；最终提供 HTTPS 链接，不仅提供本地路径。
