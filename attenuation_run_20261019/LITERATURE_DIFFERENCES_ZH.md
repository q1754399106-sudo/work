# 方法归属与本次实现范围

检索日期：2026-10-19。仅公开书目信息与方法关键词用于查询，未上传研究包内容。

| 文献 | 可核对的依据 | 本次候选与其关系 |
|---|---|---|
| Ma et al., 2020, *Multichannel absorption compensation with a data-driven structural regularization*, Geophysics 85(1), V71–V80, DOI 10.1190/geo2019-0132.1 | 原任务书给出的书目；本次 SEG 原站抓取失败，未完成原文复核 | 多道结构正则化已有先例；保留 Frozen 原有 H，不将其称为本次创新。 |
| Ma et al., 2022, *Inversion-based attenuation compensation with dip constraint*, Petroleum Science 19, 543–553 | [期刊机构站全文](https://www.cup.edu.cn/petroleumscience/docs/2022-05/5d1e01f009fe454bbbee942967d3fa9c.pdf) | 该工作在反演中引入地震倾角信息。本次 P1 不增加倾角估计，重点为固定先导下的噪声与 Q 扰动误差预算。 |
| Kaipio & Kolehmainen, 2013, *Approximate marginalization over modelling errors and uncertainties in inverse problems* | [Oxford Academic 章节摘要](https://academic.oup.com/book/12043/chapter-abstract/161414792)；搜索索引可读摘要，直接全文访问失败 | 近似误差分布已有理论基础。本次有限 bootstrap 的均值/收缩协方差是特定的固定先导代入近似，不等于精确边缘化，也不保证对地震补偿有效。 |
| Calvetti, Dunlop, Somersalo & Stuart, 2018, *Iterative updating of model error for Bayesian inversion*, Inverse Problems 34, 025008 | [作者预印本](https://arxiv.org/abs/1707.04246)、[作者机构记录](https://authors.library.caltech.edu/records/ekz96-xk362) | 其方法迭代更新模型误差分布；本次固定先导、固定协方差，一次构造后求二次问题，不借用其特定收敛结论。 |

本次工程实现采用全道集低秩 Q 误差协方差，保留跨道相关性；噪声方差上分位为保守代入，不能写作置信区间。频率仅通过 Q/子波联合校准影响 Q 扰动，不直接加入以带限波形为未知量的 G。

候选 Split 使用完整的 `M*Gref` 及传播后的协方差；历史 Frozen Split 使用近似分解算子。跨版本比较因此同时包含误差预算与路线定义变化，须结合候选自己的 Direct 和控制消融解释，不将结果归因于“分步本身更优”。

图形流程参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. [arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065)。本次核对最新记录为 2026-09-02 修订；这是制作流程来源，不是算法有效性的科学证据。
