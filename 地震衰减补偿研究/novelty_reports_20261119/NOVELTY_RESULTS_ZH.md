# 三份正式查新报告的登记与处理（2026-11-19）

三份报告由外部查新完成，对应 11-19 16:41 提交的三份查新提示词（6c75d2f）。本文件登记评级、处理办法和待一手核验的条目。报告原文为同目录下的三份 PDF。

| 文件 | 对应研究线 | 提示词 |
|---|---|---|
| `qscan_wavelet_anchoring_novelty_report.pdf` | 扫描 Q 的子波锚定理论（`qscan_anchor_theory_20261030`） | `qscan_anchor_theory_20261030/NOVELTY_PROMPT_ZH.md` |
| `avo_q_calibration_novelty_report.pdf` | 叠前 Q 补偿的 AVO 误差预算与等长窗谱比自校准（`uq_avo_20261118` + `offset_q_avo_20261118`） | `offset_q_avo_20261118/NOVELTY_PROMPT_ZH.md` |
| `q_uncertainty_novelty_report.pdf` | Q 不确定度下的期望补偿算子（叠后 `qunknown_comp_20261109`，叠前 `uq_rtm_20261117`） | `uq_rtm_20261117/NOVELTY_PROMPT_ZH.md` |

SHA-256：见 `SHA256SUMS.txt`。

评级含义（沿用报告）：A 级 = 检索范围内未见相同报道；B 级 = 部分相近；“相同”= 已有先例或已是共识。三份报告都说明检索范围不含 WoS 与 CNKI 闭库全文，部分条目只读到摘要或题录。

## 1. 扫描 Q 子波锚定理论

| 贡献 | 内容 | 评级 | 处理 |
|---|---|---|---|
| 1 | 闭式选择方程，α_a = a₀⟨a⟩/⟨a²⟩（初稿旧记号 κ），压缩斜率 1 − α_a | B | 物理本质与地震学 Shearer 等 (2019, 2022) 的震源谱–t* 权衡同源。写成“把地震学已知的权衡解析化为勘探扫描 Q 的闭式选择方程”，不写“首次发现” |
| 2 | 运行前提交数值预测，换几何、深锚窗、两工区实测逐一检验 | A | 保留 |
| 3 | 两锚点诊断：两个深度的子波窗各选一次 Q，反解 α_a 与 β | A | 核心 |
| 4 | 调稀疏权重只平移 β；迭代反馈按 α_a 收敛并放大偏差；重估子波有曲率代价 | A | 核心 |
| 5 | lasso 先验与密集反射系数失配造成无噪声下的确定性 Q 偏差 | A | 核心。程亮等 (2025, 《地球物理学报》68(2), pp. 697–698) 写明该现象“尚未出现严谨的理论证明”，可作引言切入点 |
| 6 | “残差 + 稀疏度”准则的偏移随模型而变，常数去偏不可行 | A | 核心 |

- **体裁：** 报告建议写完整长文（Geophysics 正文或《地球物理学报》），不拆成快报。
- **审稿薄弱点（本仓库自查，非报告内容）：** 实测全工区 26 条测线只在默认子波窗上测得 α̂_a = 0.681（`PAPER_DRAFT_ZH.md` 表 1）。两锚点实测只有 10 条测线、两个窗深，没有展示 α̂_a 随子波窗深度连续变化、跨过 α_a = 1 后选择反向的曲线。下一步的滑动子波窗实测（t₀ ∈ [1.4, 2.6] s）需先写冻结方案再运行。
- **待补图：** `PAPER_DRAFT_ZH.md` 图 3（加截短窗点）、图 6（加新点）、图 8（待做）。

## 2. 叠前 Q 补偿的 AVO 误差预算与等长窗自校准

| 创新点 | 内容 | 评级 | 处理 |
|---|---|---|---|
| 1 | 无拟合参数的闭式 AVO 增益 Γ(x)，以及含 B/A 的 Δ(G/A) 一阶灵敏度 | A | 核心 |
| 2 | AVO 梯度误差 ∝ σ_lnQ，定量区分算子改进与 Q 精度的贡献，给出 σ_lnQ 临界值（0.2–0.45） | A | 核心 |
| 3 | 期望补偿算子因远道 t* 绝对不确定度更大而自造负的假 AVO 梯度 | A | 核心 |
| 4 | 未动校 CMP 沿时距曲线取不等长窗造成时差压缩、谱抹平和正偏差；远道窗须与近道等长 | A | 核心 |
| 5 | 谱比自校准 + 沿路径补偿；“只按垂向 t* 补偿破坏 AVO” | B（“垂向破坏 AVO”已是共识） | 垂向补偿只作基准对照，引用 Luh (1993)、Li 等 (2015, Figs. 11–13)、程万里等 (2023)。“10.7 倍”不再作为新发现。流程与 WO2002065372A1、Geophysics 91(4):S93 部分相近 |
| 6 | 频率饱和阻尼算子 | B（广义 Tikhonov L2 的物理参数化变体） | 移出核心创新点，只作对比算子之一 |

处理已写入 `offset_q_avo_20261118/PAPER_OUTLINE_ZH.md`（主张重排）。

## 3. Q 不确定度下的期望补偿算子（叠后 + UQ-RTM）

| 创新点 | 内容 | 评级 | 处理 |
|---|---|---|---|
| 1 | 期望 Wiener 型反 Q 算子 H = S·E[G]*/(S·E\|G\|² + N) | B | 公式本身源自 Slepian (1967, JOSA 57(7), 918–922, Eqs. 12, 17–18)。写成“把随机传递函数的期望 Wiener 滤波引入含 Kolsky 频散的地震反 Q 补偿”，不写“首次提出” |
| 2 | 由 Q 不确定度与 t* 解析导出时变频带 f_Q(t)，替代人为 G_lim / f_c | A | 核心 |
| 3 | 用指数算子对随机 1/Q 的 Jensen 不等式统一解释“过补偿比欠补偿更伤” | A | 核心，但须写入 11-19 的左尾机制（见下） |
| 4 | 受限贝叶斯 + 单侧截断先验（CB/CBL） | B | 叠后线封盘，CBL 不再发表（第 4 节） |
| 5 | Q-RTM 成像条件逐点施加 Q 期望算子 + 程函 t* + 两遍随机极性成像域噪声估计 | A | 核心，但噪声估计有已知缺陷（见下） |

### 3.1 报告没有看到的 11-19 晚间结果

提示词提交于 16:41，此后两项结果改变了上表的可用程度：

1. **`uqw_rho_ablation_20261119`（20:34）：** 噪声主导（FD 14、8 dB）时，UQ-W 相对 VDB 的优势几乎全部来自对数正态先验在 u ∈ [0.1, 0.5] 的左尾（截断到 [0.5, 2] 损失 0.081 / 0.079，6/6；截断到 [0.1, 2] 恢复到 0.001 以内）；信噪比合适时（半合成）左尾只有约 0.011。所以“Q 只知道在 ×2 以内”不是机制的准确表述，Jensen 不对称的论述必须包括左尾。
2. **`uq_rtm_20261117` 阶段 4 与事后检查（20:28）：** Volve F-5 井旁相关，所有补偿方法都低于不补偿 R0 = 0.63。f_ref 改 25 Hz 后 UQW1 为 0.37–0.56，只补振幅为 0.50–0.59。UQ 在补偿方法中领先，但没有跨过“不低于不补偿”。成像域随机极性噪声估计把画弧当作信号，噪声项在实测上不起作用；阶段 3 也显示它在横向相关噪声下失效（Q 准时 −0.04 至 −0.09）。

因此创新点 5 的“两遍随机极性噪声估计”目前只适用于各道独立的随机噪声。在修好相干噪声下的成像域信噪比估计之前，不能作为无条件的贡献去投稿。

### 3.2 引文更正

| 原记录 | 更正 | 位置 |
|---|---|---|
| Zhao, Mao & Ren 2018，*Geophysics* 83(1) S15 | Wang, Y., Zhou, H., Chen, H., & Chen, Y. (2018). Adaptive stabilization for Q-compensated reverse time migration. *Geophysics*, 83(1), S15–S32. doi:10.1190/GEO2017-0244.1 | `uq_rtm_20261117/LIT_PRECHECK_ZH.md`（已改）；`NOVELTY_PROMPT_ZH.md`（加注，原文保留） |
| Eldar & Merhav 2005，*IEEE Trans. Information Theory*，doi:10.1109/TIT.2005.851730 | Eldar, Y. C., & Merhav, N. (2005). Minimax MSE-ratio estimation with signal covariance uncertainties. *IEEE Trans. Signal Processing*, 53(4), 1335–1347. doi:10.1109/TSP.2005.843701。内容是最小最大 MSE 比值后悔，不是受限贝叶斯 | `qunknown_comp_20261109/DESIGN_ZH.md`（已改）；`PLAN_MMR_ZH.md` 是运行前提交的方案，不改，以本表为准 |
| 某 CSEG 2009 摘要（作者、题名未记） | Changjun Zhang (2009)，CSEG / GeoConvention 2009 摘要；报告据此追到全文 Zhang, C., & Ulrych, T. J. (2007). Seismic absorption compensation: A least squares inverse scheme. *Geophysics*, 72(6), R109–R114。贝叶斯先验加在反射系数上，不在 Q 上。2009 摘要的题名报告未给出，待核 | `uq_rtm_20261117/LIT_PRECHECK_ZH.md`（已改） |
| （缺）期望 Wiener 滤波的出处 | Slepian, D. (1967). Linear least-squares filtering of distorted images. *J. Opt. Soc. Am.*, 57(7), 918–922. doi:10.1364/JOSA.57.000918 | `qunknown_comp_20261109/DESIGN_ZH.md` 第 2.4 节（已补） |

## 4. 叠后 UQ-W / CBL 线封盘（2026-11-19）

**决定：** 叠后 UQ-IQ / UQ-W / CBL 不单独成文，也不再开新实验（不建第七批半合成，不做“左尾 × 噪声项”交互消融）。叠后的理论部分（Slepian 期望 Wiener 引入、解析 f_Q(t)、Jensen 凸性与左尾机制）并入 UQ-RTM 论文，作为“一维理论机制”一节；CBL 舍弃。

依据：

- 新颖性：叠后期望算子（创新点 1）与 CB/CBL（创新点 4）都是 B 级，弱于叠前部分。
- 实测：`field_qperturb_20261113`，Volve 105 条线、8 组井，9 个 m 平均 UQ-W − VDB = −0.071，Q 偏大一侧 −0.135，方向与半合成相反。
- 合成确认：`uqw_shaping_confirm_20261119` 中 Q 接近正确时 UQ-W + 时不变整形比仅整形低 0.014（0/6），UQ 部分只在 Q 给错时有贡献。
- 机制：`uqw_rho_ablation_20261119` 表明噪声主导时的优势来自左尾，不来自“×2 以内”。

已有结论和文件全部保留，不删除；标签相应改为“已封盘”（`目录标签.csv`、`README.md`）。

## 5. 后续阶段（计划，未开始）

| 阶段 | 内容 |
|---|---|
| 1：子波锚定理论投稿 | 写滑动子波窗实测冻结方案（Volve 26 条线，t₀ ∈ [1.4, 2.6] s），运行；补图 3、6、8；核验第 6 节条目；定稿 |
| 2：AVO 误差预算投稿 | 核验 Geophysics 91(4):S93 与 WO2002065372A1 原文；补表 1、表 3；按目标期刊决定是否补 Volve 叠前道集“不等长窗 vs 等长窗”实测或偏移后角道集 |
| 3：UQ-RTM 攻坚 | 相干噪声下的成像域信噪比估计（新模型盲测）；UQ-RTM 先验左尾消融；在更宽频带实测上做预登记检验，事先写定 f_ref（取偏移频带）或只补振幅，目标是不低于不补偿 |

## 6. 写入参考文献前必须一手核验的条目

报告中以下条目只读到摘要、片段或题录，或题录信息不完整：

**子波锚定：**
- Cheng, L. et al. (2025). *Geophysics* 90(5), R275–R284. doi:10.1190/geo2024-0746.1（只读正文片段；需核实测子波窗位置）
- Aghamiry, H. S., & Gholami, A. (2018). *Geophysics* 83(4), V233–V242. doi:10.1190/geo2017-0001.1（片段）
- Wang, H. et al. (2021). Q-compensated denoising of seismic data. *IEEE TGRS* 59(4), 3580–3587. doi:10.1109/TGRS.2020.3010813（题录/他引）
- Ma, M. et al. (2018). *J. Appl. Geophys.* 159, 434–445（题录/他引）
- Li et al. (2026)：报告称引用 Cheng 2025 并用 Cauchy 范数防过稀疏，但没有给出题名、刊名和 DOI
- Shearer, Abercrombie & Trugman (2022), *JGR* 127(7), e2022JB024243（题录，经 Bindi 等 2025 复现）；Bindi 等 (2025, BSSA) 题录不完整
- Lari & Gholami (2019), *Geophysics* 84(1), V1–V9；Margrave (1998) SEG 摘要（片段）
- 报告注明受付费墙未核：Cheng et al. 2019、2020、2021，Gholami 2015，Chai et al. 2019

**AVO 误差预算：**
- “Obtaining effective Q factors at different offsets from CMP records”，*Geophysics* 91(4), S93（2025/2026，只读摘要/片段；作者、年份待核）
- WO2002065372A1（ExxonMobil 2002，Method for spectral balancing offset seismic data；只读摘要）
- US11294087B2 说明书、Aharchaou & Neumann (2019) *Geophysics* 84(5) V295–V305（摘要）
- Du 等 (2022) *Earth Planet. Phys.* 6(3) 259–274（摘要）
- Trantham & He (2017)、陈哲等、李国发等 (2014)（报告建议补查全文）
- Braga & Moraes (2013)、Chen et al. (2018, IEEE GRSL 15:1667–1671)（只核实到题录）

**Q 不确定度：**
- Changjun Zhang (2009) CSEG 摘要题名
- Arbel 等 (2003) 大气 Wiener 滤波；Ward & Saleh (1985)；Guan & Ward (1989)（报告只在问题 4 中提到，未列入比对表）
- Muhumuza 等 (2019) arXiv:1911.01192v1（预印本，需查正式发表版本）
- CNIPA（IPC G01V1/30）与 CNKI/万方学位论文库闭库复核

标签：**报告**（查新登记与处理决定，无新计算）。
