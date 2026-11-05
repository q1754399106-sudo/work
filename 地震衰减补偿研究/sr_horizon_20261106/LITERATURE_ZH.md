# SR 估 Q 的颜色偏差：查新记录

日期：2026-11-06。**覆盖范围有限：** 只查了公开网页（网页搜索、开放的摘要页和部分正文）。WoS 和 CNKI 没有查，用户登录学校认证有困难。大多数文献只读到摘要或检索摘要，没有读全文。所以下面的“没找到”不能算“没有”，投稿前必须在 WoS、Scopus、CNKI 补查。

## 已有工作

| 内容 | 文献 | 读到的程度 | 与本方向的关系 |
|---|---|---|---|
| 用近井反射序列的局部谱校正反射谱，再做 QVO | [Hackert & Parra 2004, Geophysics 69(6):1521–1529](https://repository.geologyscience.ru/handle/123456789/39059?show=full) | 摘要 | 测井校正的先例，方差和负 Q 减少。要有近井 |
| 用测井反演的反射序列校正频谱，非平稳 Q 模型 | [Li, Zhou, Zhao & Marfurt 2016, JSE 25:447–461](https://www.geophysical-press.com/journal/JSE/articles/233) | 检索摘要 | 同上 |
| 对数谱比公式里加反射系数校正项 | [CREWES 2012](https://www.crewes.org/Documents/ResearchReports/2012/CRR201211.pdf) | 检索摘要 | 说明白反射假设和测井校正已是常识 |
| 反射系数的幂律谱 | [Phythian 等，综述 Walden–Hosken](https://connectsci.au/eg/article/26/3/497/79934/Spectral-characterisation-of-reflectivity) | 检索摘要 | 颜色本身早已熟知；指数约 0.5–1.5，地层间有差别 |
| Gabor 反褶积中校正白反射假设 | [Montana & Margrave 2006](https://old.cseg.ca/assets/files/resources/abstracts/2006/172S0131.pdf)、[Cheng 2013](https://prism.ucalgary.ca/handle/11023/871) | 检索摘要 | 针对反褶积，不是 Q 偏差 |
| 同一反射层跨炮检距的谱比（QVO） | [Dasgupta & Clark 1998, Geophysics 63(6):2120](https://next.sciencegate.app/document/10.1190/1.1444505) | 检索摘要 | 思路最近：同一反射层抵消反射项，但用叠前道集 |
| 叠前联合反演，减少频谱干涉 | [Reine, Clark & van der Baan 2012, Geophysics 77(1):R45](https://sites.ualberta.ca/~vanderba/papers/RCV12a.pdf) | 检索摘要 | 同上 |
| 对数幅度谱的层间 Q 反演，β 为振幅标量 | [US 7,376,517（Rickett，Chevron）](https://patents.google.com/patent/US7376517B2/en) | 正文前段 | 模型和多窗 SR 同类；**β 不是颜色指数**，也没有反射谱项 |
| 联合估计反射系数和 Q（字典学习） | [Yan 等 2023](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2023.1121956/xml) | 检索摘要 | 反演类，仓库的实验显示 Q 信息会被反射系数吸收 |
| 反射数据 Q 层析 | [US 10,317,548](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10317548)、[ASEG 2015](https://publish.csiro.au/ex/ASEG2015ab264) | 检索摘要 | **可能隐含相近做法，未核对**，是新颖性的最大风险 |

[Raji & Rietbrock 2013, J. Geophys. Eng. 10(4):045012](https://academic.oup.com/jge/article/10/4/045012/5110337) 的可读部分没有把反射谱颜色作为偏差讨论。[SEP 报告](https://sepwww.stanford.edu/data/media/public/docs/sep147/yishen1/paper.pdf) 取到的 PDF 没有可读文字，这条没有核对。

## 公开渠道没找到的

1. 直接量化**随深度变化的颜色**对叠后 SR 估 Q 偏差的工作。
2. 把颜色指数和 Q 联合估计的方法。
3. 把颜色作为扰动量算 Q 区间的方法。
4. 叠后数据上比较**同一地层在不同旅行时**的频谱，以抵消任意形状的颜色。
5. 在“所有窗共用一个谱”和“每个地层一个谱”之间做部分合并（partial pooling）的 SR。

## 待补查（用户能登录时）

WoS 建议检索式：

```text
TS=(("Q estimation" OR "quality factor" OR attenuation) AND (seismic) AND (reflectivity NEAR/3 (colo$r OR spectr* OR "non-white")))
TS=(("Q" OR attenuation) AND seismic AND (horizon OR stratigraph*) AND ("spectral ratio" OR "log spectral") AND (post-stack OR poststack OR stacked))
TS=("Q tomography" AND reflection AND (reflectivity OR "spectral feature"))
```

还要读全文的：Hackert & Parra 2004 和 Reine 2012 Part 1，看它们怎样处理反射项；Q 层析的两篇，看是否已用到同层不同路径来消除反射谱。

CNKI 检索词：沿层 Q 估计、层控 品质因子、反射系数 谱 校正 Q、吸收 层析 叠后。
