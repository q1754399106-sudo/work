# 公开文献初查（不是完整查新）

日期：2026-11-06。只查了公开网络。**WoS、CNKI 没有查，所以下面的“没找到”不等于“没有”。** 标明“全文读过”的才算读过。

## 和本工作直接相关

| 文献 | 读到什么程度 | 做了什么 | 和本工作的关系 |
|---|---|---|---|
| Ratnett, Baldock, Masoomzadeh, 2019, SEG 年会摘要, doi 10.1190/segam2019-3214689.1 | **全文读过** | 用合成事件定性展示：反 Q 的增益限制会扭曲 AVO。做法：逐点 Q 场和频率上限场（取窗口峰值频率），上限内不限增益做反 Q，再按上限高切 | **最接近的先例。**“增益限制扭曲 AVO”和“限制频率而不是增益来保 AVO”的思路已经发表，**不能作为我们的新发现** |
| Li, Liu, Zheng, Huang, 2015, Absorption decomposition and compensation via a two-step scheme, Geophysics 80, V145–V155 | 只见 Du 等的引用和概述，全文未读 | 把吸收分解为时间和炮检距两部分，分两步补偿 | **分步法本身的先例**，需要读全文 |
| Removing absorption effect on AVO via inverse Q filtering in offset direction | 只见搜索摘要，ResearchGate 返回 403 | 先炮检距后深度的两步法，报告炮检距一步校正了大部分 AVO 趋势 | 分步法先例，用户需提供全文 |
| Du, Wang, Sun, Fu, 2022, Seismic attenuation compensation with spectral-shaping regularization, Earth and Planetary Physics, doi 10.26464/epp2022024 | 读了摘要、引言、方法部分 | 写出一般形式 (AᵀA + λDᵀD)，用基于数据谱的对角整形算子，目标是保低频。前向算子随炮检距变化 | **一般形式 A'A + λD'D 是已知的。** 它的算子不是按炮检距路径设计的，读到的部分也没有定量 AVO 评价 |
| Braga, Moraes, 2013, High-resolution gathers by inverse Q filtering in the wavelet domain, Geophysics | 只见搜索摘要 | 小波域反 Q 用于道集，报告 AVO 可保持 | 需要作为基线文献阅读 |

## 还没读到（403 或仅见标题）

- Prestack Q compensation with sparse tau-p operators（ResearchGate 403）：搜索摘要显示前向是频率域，用 Tikhonov，再联合去噪。
- Prestack data attenuation compensation based on inversion，EAGE 2021（earthdoc 403）。
- Chen et al., 2018, Data-driven attenuation compensation via a shaping regularization scheme, IEEE GRSL 15；Wang et al., 2015，地球物理学报 58。
- 专利 US 11294087（方向性 Q 补偿）只见摘要：远炮检距一维反滤波会过度补偿。

## 对“创新点”的影响

**不能主张的：**
- “增益或正则会扭曲 AVO”这个现象本身（Ratnett 等已有）。
- 炮检距、深度分步补偿（Li 等已有）。
- 广义 Tikhonov 这个形式（已是常规）。

**目前可以主张、但都有待 WoS 和 CNKI 核实的：**
1. **闭式偏差公式**：无拟合参数定量预测 Tikhonov 型补偿造成的假 AVO 梯度（新模型上斜率 0.986–1.027，R² ≥ 0.977）。读到的文献都只是定性展示。
2. **一步形式的频率饱和正则算子**：L 由炮检距路径的吸收构造，在 f_n 以下 ρ(f, x) = ρ(f, 0)，f_n 由闭式噪声上界决定。和 Ratnett 等的差别：这是正则化反演，不做高切，f_n 以上的信号保留；上限由噪声上界算出，且随炮检距变化，不是取峰值频率；有明确的传递函数证明。
3. **分步法与广义 Tikhonov 的频率响应对应关系**（`DIAG_ZH.md`），可以把分步法放进统一框架。
4. **对评价指标的说明**：全局归一化的 M_G、E_G 混入时变分辨率。

**需要诚实写明：** 它和分步法在机制上相近，优势主要是一步求解、波形更好、不靠抵消；对 NMO 拉伸和速度误差没有帮助。

## 来源

- [Ratnett et al. 2019](https://www.tgs.com/hubfs/Ratnett_segam2019-3214689.1.pdf)
- [Du et al. 2022](http://www.eppcgs.org/fileDQYXXWL/journal/article/dqyxxwl/newcreate/epp2022024.pdf)
- [Removing absorption effect on AVO via inverse Q filtering in offset direction（ResearchGate，403）](https://www.researchgate.net/publication/285025852_Removing_absorption_effect_on_AVO_via_inverse_Q_filtering_in_offset_direction)
- [Braga, Moraes 2013（Semantic Scholar）](https://www.semanticscholar.org/paper/High-resolution-gathers-by-inverse-Q-filtering-in-Braga-Moraes/2346c4fd3996727c69817ded472be4c87808b0c0)
- [Prestack Q compensation with sparse tau-p operators（ResearchGate，403）](https://www.researchgate.net/publication/333944283_Prestack_Q_compensation_with_sparse_tau-p_operators)
- [Prestack data attenuation compensation based on inversion（earthdoc，403）](https://www.earthdoc.org/content/papers/10.3997/2214-4609.202112741)
- [专利 US 11294087](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11294087)
