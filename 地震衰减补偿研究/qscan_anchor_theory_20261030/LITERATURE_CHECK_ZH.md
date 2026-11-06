# 文献核对：原实现能否复现压缩、震源谱与 t* 的权衡（2026-11-02）

对应投稿前四项缺口中的第 1、2 条。β 追查按建议暂停。核对方式：用 Crossref 逐篇核对作者、卷期和 DOI，摘要能读到的读了摘要；GitHub 仓库与代码搜索；能打开的开放全文读了相关段落。每条标出核对程度：**[全文]** 读过相关段落，**[摘要]** 只读摘要，**[转述]** 只见于他人转述或检索摘要，没读原文。

## 结论

1. **第 1 条目前做不了"原实现复现"。** 已发表的扫描 Q / 稀疏度选 Q 方法都没有公开代码，数据声明只写了"向通讯作者索取"。唯一公开的相关原实现是 Wang 等 (2021) 的 Q 补偿稀疏反演（GPL），但它的 Q 和 Ricker 主频都由人工给定，没有选 Q 这一步。
2. **第 2 条：地震学对震源谱与 t* 的权衡有大量正面讨论**，结构和本文相同：数据只能定出相对量，绝对水平要靠一个外加约束来定（固定小事件的角频率，或假设源谱形状）。勘探方面只有定性的说法（质心频移法、峰值频移法依赖假设的源谱形状；子波误差会使谱比法的 Q 有偏）。**没找到**把"子波窗的假设衰减"定量映射到扫描 Q 选择的工作。
3. **对论文的影响：** κ 公式不能再写成"未见讨论的新偏差机制"。它是这一经典权衡在"浅层锚定＋全局 Q 扫描"中的具体形式，数学上是初等的：反射资料只约束 t* 的深浅之差 (a(t) − a₀)u；用单一全局 Q 的 v·a(t) 去拟合 u·a(t) − a₀(u − v_a)，过原点回归就得到 κ = a₀⟨a⟩/⟨a²⟩。和规则 A + K = 1（平移不变）对应 Shearer 等 (2019) 的"相对量可分辨、绝对量不可分辨"；Shearer 等 (2022) 固定小事件角频率，对应本文的锚点。能站住的贡献要改写为：(1) 该权衡在这类勘探方法中的闭式选择方程，以及运行前提交的预测检验；(2) 两锚点诊断；(3) 说明调稀疏权重、迭代反馈、重估子波都绕不开它，因为绝对水平本来就不可辨。新颖性主要在应用和诊断上，理论本身不新。还要等第 1 条证明已发表方法确实受影响，文章才成立；比较现实的定位是短文或快报。

## 1. 已发表方法的原实现

### 1.1 方法与代码可得性

| 文献 | 选 Q 方式 | 子波 | 代码 |
|---|---|---|---|
| Cheng L., Chen M., Zhang Z., Liao J., Xu C., Li A., 2025. Q-estimation and compensation for attenuated seismic data with the scanning-Q nonstationary sparse inversion method. *Geophysics* 90(5), R275–R284. doi:10.1190/geo2024-0746.1 [摘要] | 扫描 Q，最小化浅、深补偿结果的相对误差 | 摘要只写"a suitable seismic wavelet"；SEG 页面 403，全文未读 | 未公开 |
| Aghamiry H. S., Gholami A., 2018. Interval-Q estimation and compensation: An adaptive dictionary-learning approach. *Geophysics* 83(4), V233–V242. doi:10.1190/geo2017-0001.1 [摘要] | 分组逐层搜索最优 Q，使系数最稀疏 | 未读 | 未公开 |
| Ma M., Wang S., Yuan S., Gao J., Li S., 2018. Multichannel block sparse Bayesian learning reflectivity inversion with l_p-norm criterion-based Q estimation. *J. Appl. Geophys.* 159, 434–445. doi:10.1016/j.jappgeo.2018.09.016 [题录] | 扫描 Q，l_p 稀疏度最小 | 未读 | 未公开 |
| Cheng L., Wang H., Yuan S., Zheng X., 2019. Structural interval-Q estimation by using semi-blind time-variant sparse deconvolution. *SEG Expanded Abstracts*, 729–733. doi:10.1190/segam2019-3215626.1 [摘要] | 扫描层间 Q 模型，取 l_0.1 范数最小 | 摘要称"半盲"，细节未读 | 未公开 |
| Cheng L., Wang H., Li Y., Ji Y., 2020. Multi-trace nonstationary sparse inversion with structural constraints. *Acta Geophys.* 68(3), 675–685. doi:10.1007/s11600-020-00430-3 [题录] | — | 未读 | 未公开 |
| Cheng L., Yuan S., Wang S., 2021. Joint reflectivity and structural interval-Q estimation by using nonstationary sparse inversion. *IEEE GRSL* 18(5), 776–780. doi:10.1109/LGRS.2020.2987395 [题录] | — | 未读 | 未公开 |
| Gholami A., 2015. Semi-blind nonstationary deconvolution: Joint reflectivity and Q estimation. *J. Appl. Geophys.* 117, 32–41. doi:10.1016/j.jappgeo.2015.02.030 [题录] | 多 Q 字典 | 已知子波（据 Chai 2019 转述） | 未公开 |
| Chai X., Peng R., Tang G., Chen W., Li J., 2019. *Geophys. Prospect.* 67(8), 2003–2021. doi:10.1111/1365-2478.12838 [摘要；相关段落见 THEORY_ZH.md 第 7.1 节] | 多 Q 字典＋解析反解 | 合成例用已知子波；实测在 0.3–0.6 s 统计提取 | 未公开 |
| Yan H., Yu H., Xu T., 2023. Concurrent estimation of seismic reflectivity and Q by using an optimal dictionary learning method. *Front. Earth Sci.* 11. doi:10.3389/feart.2023.1121956 [全文] | 30 个对数等距 Q，残差＋稀疏度准则取最小 | 合成例用已知 Ricker（35 Hz）；实测"从资料浅部用盲反褶积提取"（Gholami & Sacchi 2013） | 声明"向通讯作者索取" |

检索范围：GitHub 仓库搜索（Q estimation、inverse Q filtering、Q compensation 等 11 组关键词），GitHub 代码搜索（6 组），Madagascar 源码树（只有 `user/yliu/Minvqfilt.c` 反 Q 滤波、`Mmodatten1.c` 衰减正演，没有选 Q 程序），MATLAB File Exchange 与 CREWES 工具箱（CREWES 有 Q 估计工具，但需要赞助单位账号）。

### 1.2 唯一公开的相关原实现

Wang H., Huang G., Chen W., Chen Y., 2021. Q-compensated denoising of seismic data. *IEEE TGRS* 59(4), 3580–3587. doi:10.1109/TGRS.2020.3010813。代码：`github.com/chenyk1990/reproducible_research/qcomp`（GPL v2，MATLAB，`matfun/Qcdn.m`、`time_absorp_wavelet_matrix.m`、`yc_pcg.m`）。

- 正演：时变吸收子波矩阵；子波是用户给定主频 `freqdom` 的 Ricker，不从资料估计。
- Q：用户给定常数。文中实测例取 Q = 50，并比较 Q = 60–120 后建议"从 50 开始调，看哪个最利于解释"。没有选 Q 的准则。

### 1.3 我们自己做过的不能算原实现

`RESULTS_CRIT_ZH.md` 的 Cheng 2025 型、Yan 2023 型、Ma 2018 型准则，都是在本仓库 lasso 解上按论文思路重写的，只能算"同类准则受锚定"，不能回答"原方法是否受影响"。

### 1.4 可行路径（按与"原实现"的接近程度）

1. **向作者要代码**：Cheng 2025（扫描 Q，最直接）、Aghamiry & Gholami 2018、Yan 2023。需要用户本人发邮件；可以先起草。
2. **原求解器＋文献准则**：在 Wang 等 (2021) 的原求解器外加扫描壳，准则照搬 Cheng 2019（l_0.1）或 Ma 2018（l_p），Ricker 主频按其原代码的习惯定。主频若按浅层资料定，就是假设源谱的锚定：理论预测主频（锚点）一变，选出的 Q 会按 κ 跟着移动。这比自写实现更接近原方法，但审稿人仍可能说"不是原方法"。
3. **先读全文再定**：Cheng 2025、Ma 2018、Aghamiry 2018 的子波做法需要全文（SEG、Elsevier 均无开放访问），用户如有机构权限可以下载。需要确认它们的合成例是不是用真子波：若是，κ = 0，合成例本来就测不出压缩，这可以解释为什么前人没报告。Yan 2023 已确认属于这种情况。

### 1.5 2026-11-06 补充：程亮 2025 的取子波做法与路径 2 的结果

用户无法下载 SEG、Elsevier 全文，这里只用开放来源。

- 程亮等 2025《地球物理学报》68(2):696–712（doi:10.6038/cjg2023Q0281）[全文]：第一作者与 Cheng 2025 *Geophysics* 相同，*Geophysics* 那篇摘要描述的判据也是同一种（深浅层振幅谱相对误差）；但这篇的 Q 候选来自谱比法，不是扫描。合成例用已知 Ricker 作初始子波，κ = 0。实测例取浅层 0–100 ms 的主频（约 50 Hz）定 Ricker，不做衰减校正（v_a = 0），深浅层取 1–100 与 700–800 ms。按此几何 κ_U ≈ 0.09（假设吸收从记录零时算起）。
- Yan 2023 同样是合成例用真子波、实测从浅部盲反褶积取子波。所以**这两篇的合成例测不出压缩**。
- 程亮 2025 还提到，"Q 偏差降低稀疏度"只有统计验证，没有理论证明。
- 路径 2（Wang 2021 求解器＋文献做法）已做，见 `qscan_orig_solver_20261106/RESULTS_ZH.md`：C2 型判据的斜率随子波窗加深下降（W1 1.13 → W3 0.55），但已知子波斜率 1.23 超过事先上限 1.2；C1 型判据未被锚定；W1 上两判据都没有压缩。这**不是**原实现复现，κ 的数值也没有全部符合。

- 2026-11-06 后续（`qscan_c1_factor_20261106/RESULTS_ZH.md`）：求解器不是 C1 与旧结果矛盾的原因；旧的"Cheng 2025 型"是谱斜率判据，不是式 17。式 17 在 Wang 合成数据上近似 Q_sel ≈ Q_true f_w/f₀，在旧流程半合成上 κ̂ 为正且没有平衡点。

## 2. 震源谱与 t* 的权衡

### 2.1 地震学（正面讨论，多且早）

- Scherbaum F., 1990. Combined inversion for the three-dimensional Q structure and source parameters using microearthquake spectra. *JGR* 95(B8), 12423–12438. doi:10.1029/JB095iB08p12423 [摘要]。单台窄带谱估 Q 和角频率是"strongly nonunique"：低 Q／高 f_c 与高 Q／低 f_c 拟合得一样好；用多台层析把 Q 结构作为约束。
- Abercrombie R. E., 1995. *JGR* 100(B12), 24015–24036. doi:10.1029/95JB02397 [摘要]。地表资料的"最小震源尺寸"可能是近地表强衰减造成的假象；2.5 km 井下记录分别在固定 n 和固定 Q 的模型下估参数。
- Anderson J. G., Hough S. E., 1984. *BSSA* 74(5), 1969–1993. doi:10.1785/BSSA0740051969 [转述]。高频衰减参数 κ；ω² 源模型下 κ 与 t* 等价（据检索摘要，原文未读）。注意符号冲突：论文中我们的 κ 需改名或注明。
- Ko Y.-T., Kuo B.-Y., Hung S.-H., 2012. Robust determination of earthquake source parameters and mantle attenuation. *JGR* 117(B4). doi:10.1029/2011JB008759 [摘要]。为分开源的频率依赖与介质衰减，令同簇事件共享 Q、每个事件只有一个 f_c，以减少自由度。
- Trugman D. T., Shearer P. M., 2017. *JGR* 122(4), 2890–2910. doi:10.1002/2017JB013971 [摘要]。"inherent tradeoffs in separating source and path effects"。
- Shearer P. M., Abercrombie R. E., Trugman D. T., Wang W., 2019. Comparing EGF methods for estimating corner frequency and stress drop from P wave spectra. *JGR* 124(4), 3966–3986. doi:10.1029/2018JB016957 [摘要]。全局 EGF 拟合在绝对应力降、标度和高频衰减率之间有权衡；事件间的相对谱形可分辨，绝对水平不可分辨；约束小事件 f_c 后两种方法一致。**与本文结构最接近。**
- Shearer P. M., Abercrombie R. E., Trugman D. T., 2022. *JGR* 127(7). doi:10.1029/2022JB024243 [摘要]。把附近最小事件的平均 f_c 固定为常数来消除权衡，代价是这些事件的真实空间变化也被去掉。这相当于本文的锚点。
- Shearer P. M., Prieto G. A., Hauksson E., 2006. *JGR* 111(B6). doi:10.1029/2005JB003979 [摘要]。谱分解分离源、台站、路径项，用经验格林函数校正衰减。

### 2.2 勘探

- Quan Y., Harris J. M., 1997. *Geophysics* 62(3), 895–905. doi:10.1190/1.1444197 [摘要]。质心频移法；推导需假设源谱为高斯等形状。
- Zhang C., Ulrych T. J., 2002. *Geophysics* 67(5), 1542–1547. doi:10.1190/1.1512799 [摘要]。峰值频移法，假设源谱为 Ricker。源谱形状不符时会产生系统误差，这是后续改进论文的共同出发点 [转述]。
- Rickett J., 2007. Estimating attenuation and the relative information content of amplitude and phase spectra. *Geophysics* 72(1), R19–R27. doi:10.1190/1.2399451 [摘要]。用两深度级的相位差估衰减（源谱相消），比较振幅与相位的信息量。
- Raikes S. A., White R. E., 1984. *Geophys. Prospect.* 32(5), 892–919. doi:10.1111/j.1365-2478.1984.tb00745.x [摘要]。反射时窗谱比按时窗内反射系数谱校正后更平滑。**对应 β 的反射系数谱色项。**
- White R. E., 1992. The accuracy of estimating Q from seismic data. *Geophysics* 57(11), 1508–1511. doi:10.1190/1.1443218 [题录]。
- Wang Y., Morozov I., Marfurt K., 2020. *SEG Expanded Abstracts*, 2789–2793 [转述，已在 THEORY_ZH.md 第 7 节]。谱比法依赖无衰减子波谱形的假设，子波估计误差会使结果有偏。

### 2.3 没找到的

勘探反射资料中，"子波窗假设的衰减"误差对扫描 Q 选择的闭式或定量分析（中英文都检索过）。这是本文仍可占据的位置，但要明确写成"已知权衡在这类方法中的具体形式"。

## 3. 对初稿的具体修改（2026-11-06 已并入第二版初稿；κ 符号暂保留并注明）

- 引言第三点"子波估计时假设的衰减对选出 Q 的影响，目前未见定量分析"保留，但前面要先交代 2.1 的地震学权衡和 2.2 的源谱假设问题，把 κ 定位为它的具体形式。
- 第 3.2 节补一句：κ 等于把带截距的 t* 差拟合成过原点的比例模型时的截距泄漏。
- 符号：κ 与地震学的 kappa（Anderson & Hough 1984）冲突，改名（例如 η）或注明。
- 参考文献按上表补全卷期和 DOI；Cheng 2025 的作者为 Liang Cheng, Mengxin Chen, Zhijun Zhang, Jianping Liao, Cuiyun Xu, Anyu Li。
- 2026-11-06："适用所有浅层锚定方法"不能再写，也不能把 Cheng 2025 型判据列入受影响方法，理由见 1.5。
