# UQW-S 查新记录（2026-11-15，公开渠道初查）

**范围**：网页检索，中英文共 16 组关键词。大多只读到摘要或书目信息，没有读全文。

**没有查**：WoS、CNKI、SEG Library 和 EarthDoc 的全文检索，Wang (2008)《Seismic Inverse Q Filtering》一书。投稿前必须补上。

## 逐条结论

UQW-S 有四个可能主张的新意，下表逐条对照最接近的先例。

| 主张 | 最接近的先例 | 判断 |
|---|---|---|
| 1. 对 Q（相对衰减比例 u）的不确定度取期望，得到反 Q 算子 H = E[G*] / (E\|G\|² + ρ) | **Sternad & Ahlén (1993)**：Kybernetika 29(5) 439–454 的"谨慎维纳滤波器"（cautious Wiener filter），对模型误差的概率描述取平均，设计不比普通维纳滤波器复杂；Automatica 29 (1993) 661–679；Öhrn 等 IEEE TAC 40 的多变量推广。鲁棒维纳的相关理论：Poor (1980) 最小最大，Eldar (2005) IEEE TIT 最小最大后悔 | **一般原理是已知的**，不能主张"期望维纳"本身。没有找到把它用于反 Q 的先例：带 Kolsky 频散、1/Q 取对数正态、只在子波锚点以下不确定。可主张的是**应用与推导**：σ = 0 时退化为稳定化反 Q；保守程度随深度和频率自动变化 |
| 2. 稳定项由数据的时频信噪比决定 | **Zhao, Mao & Xu (2019)**：广义稳定反 Q 滤波，J. Appl. Geophys.，按信噪比选择补偿强度和带宽；同组的时频信噪比估计。**郭志伟等 (2019)**：石油科学通报，时空变自适应反 Q，增益随局部信噪比变化。自适应增益限 (2015, 地球物理学报)；变稳定因子 (2018, JGE)；反演类方法用信噪比定稳定因子 | **已有大量先例**，单独不能主张。没有找到把 Q 不确定度与数据信噪比放进同一个维纳式的先例 |
| 3. 反 Q 后在测得的子波频带内做部分（0.5 次方）时变谱整形 | **van der Baan (2012)**：Geophysics 77(4) V133，相位反 Q 加零相位时变维纳，即本研究的 VDB 基线；时变谱白化 (TVSW)；专利中用 exp(βf)（β 约 0–0.5）的部分谱平衡 | **部分整形、时变整形都有先例。** "只在目标子波测得到的频带内整形"和"整形强度取 0.5 次方"没有找到原样先例，但属于工程处理，新意弱 |
| 4. 把对 Q 误差稳健作为设计目标 | 不需要已知 Q 的 Q 补偿稀疏反褶积（Geophysical Prospecting 2019，"Some remarks on Q-compensated sparse deconvolution without knowing the quality factor Q"）；多 Q 训练的 CNN 补偿（2024）；先稳健估 Q 再反 Q 的流程 | 目标相同，手段不同：字典、学习、估 Q，都没有对不确定度取期望。没有找到对 Q 误差做期望或最小最大设计的反 Q。**可以主张，但要和这几类方法比较或讨论** |

另外，"锚点"做法（相对浅层参考窗的衰减）在谱比法估 Q 中很常见，不单独主张。

## 对论文表述的建议

- **不要写**"首次提出对模型不确定度取期望的维纳滤波"。要引用 Sternad & Ahlén (1993) 和 Eldar (2005)，说明本文是把谨慎维纳滤波用到反 Q 补偿上。
- **新意可以写为**：
  1. 只对子波锚点以下的相对衰减取对数正态期望，得到一族反 Q 算子：σ = 0 时为稳定化反 Q，保守程度由声明的 Q 范围算出，不调参；
  2. Q 不确定度与数据信噪比在同一个维纳式中结合，消融实验显示两者缺一不可（W0、UQ-IQ 都更差）；
  3. 在目标子波测得到的频带内做部分整形，它修复了欠补偿一侧，并在两种生成器和实测井上确认、复现。
- **必须讨论的对照**：
  - Zhao 等 (2019) 的广义稳定反 Q（按信噪比选强度和带宽）。仓库里还没有实现，应作为基线补做；
  - 郭志伟等 (2019) 的时空变自适应反 Q；
  - van der Baan (2012)，即已有的 VDB；
  - 不需要已知 Q 的稀疏反褶积，与仓库的 SC、QScan＋SC 相近，但还没有与 UQW-S 比较。
- **附带的方法学发现**：VDB 基线的目标子波在噪声主导时塌缩。这与 TVSW、谱平衡文献里"过度或不足白化"的讨论有关，可以作为对比实验的公平性说明。

## 检索记录（主要来源）

- van der Baan 2012：https://sites.ualberta.ca/~vanderba/papers/VanB12.pdf
- Sternad & Ahlén 1993：https://www.kybernetika.cz/content/1993/5/439/paper.pdf；https://eudml.org/doc/27687
- Öhrn、Ahlén、Sternad 多变量谨慎维纳：https://theses.eurasip.org/theses/70/design-of-multivariable-cautious-discrete-time
- Eldar 2005：https://weizmann.ac.il/math/yonina/sites/math.yonina/files/Robust%20Deconvolution%20of%20Deterministic%20and.pdf
- Vastola & Poor：https://portal.mardi4nfdi.de/wiki/Publication:790784
- Zhao、Mao、Xu 2019（DOI 10.1016/j.jappgeo.2019.07.007）：https://www.sciencedirect.com/science/article/abs/pii/S0926985118310796
- 时频信噪比估计：https://bohrium.dp.tech/paper/arxiv/5488f23c45ce471f90964a34
- 郭志伟等 2019：https://www.cup.edu.cn/pub/sykxtb/docs/2019-06/20190628083645292032.pdf
- 自适应增益限反 Q 2015：http://www.geophy.cn/article/doi/10.6038/cjg20150727
- 变稳定因子 2018：https://academic.oup.com/jge/article/15/4/1269/5204258
- 两种抗噪反 Q 的比较：http://en.dzkx.org/article/doi/10.6038/pg2020DD0072
- 不需要已知 Q 的稀疏反褶积：https://www.earthdoc.org/content/journals/10.1111/1365-2478.12838
- 多 Q 训练的 CNN 补偿：https://www.researchgate.net/publication/383192472
- 部分谱平衡专利：https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10324211
- TVSW：https://wiki.seg.org/wiki/Time-variant_spectral_whitening_-_book
