from pathlib import Path
import json
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[1]
s=json.loads((R/'analysis/summary.json').read_text(encoding='utf-8'))
cmp=pd.read_csv(R/'analysis/predictor_comparison.csv');cv=pd.read_csv(R/'analysis/leave_model_out.csv')
cal=pd.read_csv(R/'analysis/calibration_fd_convergence.csv')
q=pd.read_csv(R/'analysis/qscan_diagnosis.csv');d=pd.read_csv(R/'analysis/derivative_checks.csv')
names={'depth_time':'双程时间（深度代理）','inverse_observed_snr':'观测窗 SNR 倒数','q_spread':'Q 扰动输出离散度',
       'U_r':'仅反射系数','U_Q':'反射系数 + 各层 Q','U_Q_f_phase':'再加入主频、相位','U_cal':'再加入独立校准（完整指标）','U_expanded':'完整指标 + 扩展支持'}
rows=[]
for _,x in cmp[cmp.subset=='live'].iterrows():
    rows.append(f"| {names[x.feature]} | {x.spearman:.3f} | {100*x.top20_energy:.1f}% |")
cvnames={'simple':'三个简单指标','simple_plus_r':'简单指标 + 仅反射系数风险','simple_plus_Q':'简单指标 + Q 风险',
         'simple_plus_f_phase':'简单指标 + Q/主频/相位风险','simple_plus_theory':'简单指标 + 完整指标',
         'simple_plus_expanded':'简单指标 + 扩展支持风险'}
cvrows=[]
for name in ['simple','simple_plus_r','simple_plus_Q','simple_plus_f_phase','simple_plus_theory','simple_plus_expanded']:
    x=cv[cv.variant==name][['spearman','r2']].mean()
    cvrows.append(f"| {cvnames[name]} | {x.spearman:.3f} | {x.r2:.3f} |")
qrows=[]
for _,x in q.iterrows():
    qrows.append(f"| {x.model}/{x.condition} | {x.selected_Q_factor:.3f} | {100*(x.q0_error_ratio-1):+.2f}% | {x.predicted_log_risk_ratio:+.3f} | {x.degradation_auc:.3f} |")
dr=[]
for name,g in d[d.step==1e-4].groupby('parameter',sort=False):
    dr.append(f'| {name} | {g.observation_relative.max():.3e} | {g.target_relative.max():.3e} |')
txt=f'''# 新推导能否提前预测失真：首轮诊断结果

**结论：部分有效，但完整新指标未通过继续开发补偿算法的门槛。** 仅反射系数的局部输出敏感性具有额外预测价值；加入 Q、主频、相位和独立校准后没有增强这一优势，且漏报了关键 QScan 退化。下一步应验证缺失的系统偏差和模型因素，暂不开发新补偿算法或增加正则项。

## Material Passport

- 类型：已有合成实验上的数值与预测诊断；已实际运行 MATLAB 和 Python。
- 模式：run；项目根目录 `D:\\ChatGPT\\研究`；日期 2026-10-24。
- 数据：primary、complex、F1、F2 四个地质模型，各 Moderate/Strong；共 8 个场景。
- 状态：导数核验通过，旧算法/旧评分复现通过，完整指标的额外预测能力未通过。历史数据已查看，不能称为未见模型的最终验证。

## 固定目标与隔离

旧 FullCalRefit、QScan 及评分未改。本轮冻结版本为 `bae_q_20261024/results`；它与更早 `complex_fieldlike_20261023` 报告的 FullCalRefit 数值略有不同，不把跨归档版本差异算成本轮改进。complex 仍以保留真实源相位、接收增益、静校正的无吸收波形为目标；F1/F2 仍以原细层模型的无吸收弹性响应为目标。**本轮没有改为统一零相位目标，也没有生成新数据。** 不同模型保留原有目标定义。

预测程序只读取 [inputs/](inputs/) 白名单文件，其中包括观测、旧 sigma、公开配置、旧输出和真实独立校准观测，不含真值。真实反射系数/源参数不用于中心解。风险文件先写出并 [SHA-256 封存](predictions/SEALED_SHA256.json)，评分进程随后读取原始真值。封存后未调整预测定义、排序或阈值。后续分析增加了消融对照及有效域标记，没有修改封存预测。输入导出是结构隔离，不能把它宣传为不知道历史结论的盲测。

每场景固定 35 道（1,4,…,97，加 23 和 78 道），0.30–1.60 s 划分为 26 个不重叠 50 ms 窗，共 **7280 窗**。主排序分析使用观测非零的 **7176 窗**；正常道定义来自观测质量，未由真值筛选。原有全道和旧正常道误差均保留在 [case_scores.csv](analysis/case_scores.csv)。窗口不是独立模型，两个噪声条件也不是两个独立地质模型；本轮不以大量窗口制造显著性结论。

## 1. 导数与数值核验

使用 z=(r,log Q₁…₄,log f,φ)，观测 A=G_Q W_f,φ，目标 B=W_f,φ。Q 对目标的导数为零，主频和相位对目标的导数非零。反射幅度包含有效接收增益，时间支持容纳时移；没有单独识别接收增益或静校正。局部源相位中心沿用旧方法的零相位假设，未读取真实相位。这里核验的是该局部工作模型，不能据此证明它包含复杂模型的真实传播与源波形。

| 参数/方向 | 观测导数最大相对误差 | 目标导数最大相对误差 |
|---|---:|---:|
{chr(10).join(dr)}

表取中心差分步长 10⁻⁴，最大值跨所有场景。步长还检查了 10⁻²、10⁻³、10⁻⁵；Q 的零目标导数绝对误差为零。观测/目标的最坏相对误差分别为 **{s['derivative_observation_max_at_1e4']:.3e} / {s['derivative_target_max_at_1e4']:.3e}**，低于预定 10⁻⁵ 标准。

独立校准使用原 6 个振幅/正交相位模板，先剖除其自由度，再由校准观测直接形成 K，不使用旧协方差代理。只有 complex 有此独立观测；其余模型明确无 K。校准对每事件相位自由，因此不能约束剖面局部源相位；局部 K 数值秩为 5。白化拟合残差 RMS 为 1.527，提示工作校准模型也有未解释变化；不能把这项信息直接称为已校准的概率保证。

校准 K 本身采用中心有限差分；补充检查 10⁻⁴ 相对 10⁻⁵ 的导数收敛，代表性基准/Moderate QScan 中心最大差 **{cal.loc[cal.step==1e-4,'relative_to_step_1e5'].max():.3e}**。这属于数值步长稳定性检查，不冒称独立解析推导核验。

原始算法重算相对归档最大差：FullCalRefit **{s['frozen_base_max']:.3e}**，QScan **{s['frozen_scan_max']:.3e}**，与保存为 single 的精度一致；旧正常道评分最大绝对差 **{s['old_score_max_difference']:.3e}**。新诊断零相位正演与旧未归一化字典一致。数值 SVD 截断 10⁻⁸、10⁻¹⁰、10⁻¹² 在有效观测窗给出相同风险；未发现数值零空间导致的无穷输出风险。

![导数收敛](figures/derivative_convergence.png)

## 2. 能否提前排序实际误差

风险为点态 sqrt(g H† gᵀ) 的窗内 RMS，采用旧 sigma 白化、epsilon_loc=1。它用于比较局部噪声/参数敏感性，**不是 95% 区间，也不是已覆盖全部失配的真实误差上界**。对照 Q 离散度来自同一旧算法在 Q±10%、±30%、±50%、Q0 的输出标准差。深度使用双程时间代理；SNR 仅用观测窗 RMS/旧 sigma，未用真实噪声或真实信号。

| 预测指标 | 场景平均 Spearman | 风险最高 20% 窗捕获的实际误差能量 |
|---|---:|---:|
{chr(10).join(rows)}

8 场景等权平均，均为正常道。逆 SNR 与绝对误差的负相关不等于它没有信息：低能量窗的绝对波形误差也较小。下节回归允许各简单指标自由取正/负系数，避免靠选错 SNR 方向制造新指标优势。

完整指标的相关为 0.546、能量捕获 34.6%，略高于 Q 扰动离散度的 0.527、34.2%，但远低于仅反射系数风险的 0.639、41.4%。**原始相关略高并不足以证明附加预测能力。** 加入干扰参数使允许变化变大，并不保证其排序更贴近实际算法误差。

![风险与实际误差的直接对应](figures/risk_error_correspondence.png)

灰点是各已有场景的风险十分位实际 MSE，纵轴除以该场景平均 MSE；蓝线为 8 场景等权均值。它直接显示高风险窗是否更易出错，未作为调参依据。底表：[risk_error_deciles.csv](analysis/risk_error_deciles.csv)。

![预测指标对照](figures/predictor_comparison.png)

![空间对应](figures/risk_error_maps.png)

complex 的第 23、78 道是零观测，旧 sigma=realmin。此时白化及固定支持风险不满足有效观测条件，封存计算中的零值是无效占位，**不能解释为可靠或有限误差保证**；扩展候选也不能修复该噪声/观测模型。图中风险位置标灰，实际误差保留。原值未删除或改写，有效域另见 [prediction_domain.csv](analysis/prediction_domain.csv)。全道旧评分照常计入其损失。

## 3. 是否超出简单指标的预测能力

用固定线性秩回归：基础输入为时间、观测 SNR 倒数、Q 扰动离散度；再分别加入各风险消融。每次同时留出一个地质模型的两种噪声，输入只在各场景内转分位秩。没有随机打散同一剖面的窗口，没有搜索回归超参数。以下为 4 个留出地质的等权平均。

| 输入 | 留出模型排序相关 | 留出模型秩回归 R² |
|---|---:|---:|
{chr(10).join(cvrows)}

仅反射系数风险将相关从 0.583 提高到 0.650、R² 从 0.346 提高到 0.421，说明输出敏感性这条诊断路线有可保留的信号；但改善并不对每个地质模型一致。加入 Q/源参数后优势减弱，完整指标的相关回到 0.583、R² 降至 0.339；完整指标不提供本轮要求的额外预测证据。扩展支持也未改善平均结果。校准仅在 complex 可得，不能据此概括所有校准设计无用。

![留一地质模型预测](figures/held_model_prediction.png)

## 4. QScan 的 Q0 退化检验

| 场景 | QScan 选择的 Q 倍率 | 旧评分误差变化 | 预测 log(Uscan/Ubase) | 窗口退化 AUROC |
|---|---:|---:|---:|---:|
{chr(10).join(qrows)}

正的实际变化表示退化，正的预测分数表示预测风险增加。全部窗的风险比在有效观测域内均可计算。以分数正负预测全场景退化，**8 例仅 4 例方向正确；4 个实际退化案例只识别出 complex/Strong 1 例**。这不是事后选择最优阈值所得。

- **给定真 Q 的 primary/Moderate**：误差增加 0.70%，指标反而预测风险下降，未提前识别；其窗口 AUROC 0.619 不能替代全场景退化方向失败的事实。
- **F1/Moderate**：误差增加 16.91%，指标预测 log 风险比 -0.257，即风险下降；窗口 AUROC 0.471，漏掉关键失败。F1 的 Q0 是上尺度宏观 Q，不能称为细层全真 Q。
- complex/Strong 退化 8.29% 被识别，但 complex/Moderate 实际改善 9.42% 被误报。complex Q0 是标定 Q，不是横向真 Q。

![QScan 退化预测](figures/qscan_failure_prediction.png)

窗口 AUROC 的 8 场景平均为 {q.degradation_auc.mean():.3f}（Q 离散度对照 {q.qspread_degradation_auc.mean():.3f}）；个别案例有效，整体不足以代替关键退化预测要求。primary 的原弱反射评分也照算：Moderate 弱反射误差 0.32933→0.32684，尽管整体误差略恶化，仍说明两种评价不能混为一个“改善”。F1/F2 无对应旧弱反射定义，不补造标签。

## 5. 下一步决定

**保留：** 正确的目标导数、输出敏感性诊断及仅反射系数基线。它们确实能在读取本轮真值前产生有用排序。

**暂不进入：** 根据完整 U 调补偿强度、删模式或设计恢复频带。当前完整指标未提供额外预测能力，关键 QScan 失效也未能提前定位。

**优先补的缺口：** (1) 局部信息矩阵只限制中心附近的可变化程度，没有检验中心输出相对真实目标的系统偏移；(2) 稀疏支持与源波形中心条件可能遗漏误差方向，12 ms 固定网格扩展不等于完整多支持覆盖；(3) complex 的混合相位、延迟源分量、横向 Q 与校准迁移，以及 F1/F2 的细层传播/多次波，不由当前四层工作字典充分表达；(4) 旧标量 sigma 无法覆盖相关噪声、额外干扰与失效观测。上述是与观测及代码一致的待验证解释，本轮没有分解各自贡献。

下一项小诊断宜固定这些指标，在已有数据上逐项加入可观测的中心偏差/模型失配检验，先要求能同时区分 F1/Moderate 的退化与 F2/Moderate 的改善，再讨论算法。**不把本轮用于调试的模型重新命名为最终未见测试。** 新模型仍留待方法确定后生成。

## 复现与完整交付

- 冻结方案：[PLAN_ZH.md](PLAN_ZH.md)；理论原文快照：[THEORY_SOURCE_ZH.md](THEORY_SOURCE_ZH.md)。
- 代码：[scripts/](scripts/)；原代码与仓库版本哈希：[UPSTREAM_HASHES.json](UPSTREAM_HASHES.json)。
- 预测：[risk.csv](predictions/risk.csv)，封存：[SEALED_SHA256.json](predictions/SEALED_SHA256.json)。
- 导数：[derivative_checks.csv](analysis/derivative_checks.csv)，校准步长：[calibration_fd_convergence.csv](analysis/calibration_fd_convergence.csv)。
- 对照表：[predictor_comparison.csv](analysis/predictor_comparison.csv)，分场景：[predictor_by_scene.csv](analysis/predictor_by_scene.csv)。
- 留模型检验：[leave_model_out.csv](analysis/leave_model_out.csv)；QScan：[qscan_diagnosis.csv](analysis/qscan_diagnosis.csv)。
- 原始逐窗评分：[scored_windows.csv](analysis/scored_windows.csv)；6 组图均提供 PNG、PDF、说明文本。
- 原始执行日志：[logs/](logs/)。首次预测启动因单线程启动选项与 4 线程设置冲突而失败，日志保留为 predict_attempt1；修正启动配置后实际计算成功，未隐藏失败或伪造重跑记录。

MATLAB R2024b，Python 3.12，NumPy 2.5.4，pandas 3.0.6，Matplotlib 3.11.2。先恢复原始 `bae_q_20261024/source` 和本目录，再依次执行 `export_inputs`、`predict`、`score_outputs`；评分前由 `scripts/run_stage.py predict` 自动封存。随后运行 `scripts/check_calibration.m`、`scripts/analyze.py`、`scripts/write_report.py`。重新预测前应移走已有输出至另一个运行目录，以保留本次封存；不要覆盖后再声称是原始本轮结果。
'''
(R/'RESULTS_ZH.md').write_text(txt,encoding='utf-8')
print('Wrote RESULTS_ZH.md')
