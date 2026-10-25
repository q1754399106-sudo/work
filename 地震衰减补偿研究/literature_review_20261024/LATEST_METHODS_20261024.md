# 最新文献与方案提取（2026-10-24）

研究范围：地震非平稳反褶积、Q 补偿、子波估计与 Q 可靠性。承接本项目 bae_q 和 coherence_soft_v2 实验。本轮为定向叙事检索，不是系统综述；“最新”指此次可核实检索结果，不保证穷尽。

## 结论

优先开展“谱约束子波校正 → Q 候选选择去偏 → Q 稳定性诊断”。不再将固定 BAE 协方差或进一步调整软收缩系数作为主线。依据是本地实验报告：BAE-Q 未通过；QScan 在独立正演上稳健但因子偏高；复杂模型主要误差来自子波和正演失配。这里的机制解释仍是待检验假设，不能由现有误差曲线直接证明。

## 文献证据表

### 1. Zhao, H., Gao, J., & Chen, H. (2025)
Nonstationary blind deconvolution using spectral constraints. Geophysics, 90(4), V325–V337. DOI: 10.1190/geo2024-0284.1。

原文方法：分段平稳化，估计时变子波功率谱，再用谱约束进行盲反褶积。适合解决源子波与 Q 事先未知的问题。

核验：作者大学机构库摘要与作者代码仓库；未取得论文全文、未运行代码。机构库：https://ir.specialsci.cn/xajt/scholar/result?scholar=9a15e434ecd70556%3AJinghuai+GAO 。作者代码：https://github.com/justZHQ/Nonstationary-blind-deconvolution-using-spectral-constraints （已确认 program 目录存在，GPL-3.0；未审查具体实现）。

可提取：把固定 Ricker 字典改为由局部谱约束的时变子波字典。功率谱不能唯一确定相位，因此本项目还需显式处理相位与尺度歧义。优先级最高。

### 2. Yin, Y., Wang, C., & Wang, B. (2026)
Accurate and robust Q estimation based on low-rank constraint of the attenuation function and DBSCAN clustering. Journal of Geophysics and Engineering, 23(5), 1491–1501. DOI: 10.1093/jge/gxag053。在线发表 2026-04-09，卷期为 2026-10。

原文方法：利用衰减函数低秩性进行反演，线性拟合等效 Q；对实现参数变化产生的 Q 估计做 DBSCAN 聚类，再约束转换为分层 Q。验证包括合成、模拟 VSP 和实测 VSP。

核验：出版社摘要；搜索返回了部分表格，但未精读全文。https://academic.oup.com/jge/advance-article/doi/10.1093/jge/gxag053/8650496

可提取：多窗口、多频带、多邻道重复估计，记录主簇比例、离散度和离群率。将其用于本项目 QScan 可信度是迁移建议；聚类集中不保证物理 Q 正确，也不是概率置信区间。地面反射数据受干涉影响，不能直接套用 VSP 效果。

### 3. Geng, B., Wang, N., & Shi, Y. (2026)
Flexible Automatic Differentiation-Based Framework for Seismic Data Q Compensation. Geophysical Prospecting, 74(3), e70161. DOI: 10.1111/1365-2478.70161。在线发表 2026-03-19。

原文方法：用自动微分对正演计算图求梯度，模块化组织补偿目标与正则项。出版社摘要报告合成和实测应用；没有据此认定其解决了 Q 与子波联合可识别性。

核验：出版社摘要和数据声明；全文入口转回摘要。数据声明为联系作者合理请求获取。https://onlinelibrary.wiley.com/doi/abs/10.1111/1365-2478.70161

可提取：联合求解器的工程组织方式。当前 MATLAB 小维度子波参数可先用交替优化实现；只有参数规模扩大后才需要迁移自动微分框架。

### 4. Khaloui, F., & Haghshenas Lari, H. (2026)
Semi-blind multichannel nonstationary seismic deconvolution utilizing a skip-connected autoencoder neural network. Journal of Applied Geophysics, 250, 106192. DOI: 10.1016/j.jappgeo.2026.106192。卷期 2026-07。

原文方法：多道半盲非平稳反褶积，以粗略子波估计生成训练数据，利用跳连自编码器处理子波频率、相位、振幅变化。作者报告简单结构训练数据可用于较复杂结构例子。

核验：出版社搜索索引中的摘要和 highlights；直接页面访问 403，未读全文。https://www.sciencedirect.com/science/article/pii/S092698512600100X

可提取：子波扰动应覆盖相位、频率和振幅，而不只是 Q 扰动；神经网络作为后续对照。论文摘要不支持“任意工区泛化”，也不能推断其已解决鬼波。

### 5. Wu, J., Li, W., Ma, Y., & Li, Y. (2026)
Resolution Enhancement of Seismic Imaging via Absorption Compensation with Curvelet Transform in Frequency Domain. Applied Sciences, 16(20), 10040. DOI: 10.3390/app162010040。发表 2026-10-10。

原文方法：曲波域区分有效信号与噪声，对有效方向进行补偿增强。

核验：出版社索引的摘要与部分正文；直接打开限流，未完整阅读。https://www.mdpi.com/2076-3417/16/20/10040

可提取：作为结构域补偿对照，重点检查断层、尖灭和弱反射。它是本轮最新可核实论文，但与当前“子波失配”瓶颈的直接相关性低于前两篇；方向性增强也不能证明保幅。

### 6. Peng, J., Wang, X., Liu, Y., Li, Y., & Wang, M. (2025)
Physics-Informed Cross-Learning for Seismic Acoustic Impedance Inversion and Wavelet Extraction. arXiv:2512.11521，提交 2025-12-12。

原文方法：联合声阻抗反演与子波提取，通过物理约束交叉学习处理少井资料情形。

核验：读取作者预印本 HTML 的摘要和方法框架，未逐式审读全文。页面标注期刊名，但本轮未核实正式出版 DOI，按预印本引用。https://arxiv.org/html/2512.11521v1

可提取：子波作为待估变量、重建闭环用于约束。其任务是声阻抗反演，并非 Q 补偿；迁移到本项目不意味着方法已被验证。

## 三项可实施方案（本项目设计，不是上述论文的原样复现）

### P1：谱约束低维子波校正

输入：现有观测道集、初始子波、QScan 候选范围、噪声估计。

1. 先固定 Q，估计局部子波的主频修正与常相位；邻道采用共享或平滑参数，断层处允许分段。
2. 再以校正子波更新反射系数；最后扫描 Q。分开更新可帮助诊断 Q 是否在吸收子波误差。
3. 对子波施加单位能量、固定时间中心和极性约定，避免反射系数与子波相互任意缩放或平移。绝对振幅仍需要外部标定。
4. 第一版只开放频率和相位。鬼波时延/强度作为第二版显式项，防止一次引入过多不可识别参数。
5. 功率谱约束来自观测的稳健局部统计；不能把观测谱直接当作子波谱，因为其中还含反射系数谱与噪声。

建议优化结构：观测残差 + 反射系数稀疏项 + 子波谱约束 + 子波横向平滑项。保留 coherence_soft_v1 作为固定输出阶段，先隔离子波校正的收益。谱约束、相位参数和鬼波项分别消融。

### P2：QScan 选择去偏

动机：本地独立正演中最优因子约偏高 28%；这是本地报告观察值，不是普遍结论。

1. 保留原 QScan 作为对照。
2. 对每个 Q 候选先做稀疏支撑选择，再在相同正则规则下去偏重拟合振幅。
3. 使用连续频带留出、带缓冲的时间块留出评价预测残差；分块定义、噪声白化和频率相关性处理必须在开发阶段固定。不能把相邻样点随机拆开制造信息泄漏。
4. 原目标、去偏训练残差、留出预测残差分别比较；仅训练残差容易偏向过拟合，不能作为当然正确的新准则。
5. 若近最优候选范围很宽，输出可接受 Q 集合；不要强行声称唯一物理 Q。平均输出能否改善波形要另测，不能把“平均 Q”与“平均补偿结果”等同。

它与 P1 可分别测试：先确认偏差来自评分还是子波失配，再做组合。

### P3：Q 稳定性诊断与有条件补偿

1. 对窗口位置、窗口长度、有效频带和邻道块做扰动，重复估计 log(Q) 或 Q 校正因子。
2. 汇总主簇占比、簇内 MAD、无簇比例、多簇数量和留出残差；DBSCAN 的尺度与最小簇参数也需开发集确定。
3. 高稳定性且预测检验通过：使用 P1/P2 所选结果。低稳定性：同时给出较保守补偿及候选离散度，不输出单一可信物理 Q。
4. 以“是否能预测过补偿风险”校验诊断价值，而不是只展示聚类图。系统性模型误差可能让所有候选同时偏向同一错误解。

## 实验顺序与判断标准

先做 P1 和 P2 的独立开发消融，再组合；P3 在已有候选上增加诊断。无需先训练大型网络。

固定基线：FullCalRefit、QScanRefit、coherence_soft_v1、Wang 2006；增加 NBDSC 原作者实现作为文献对照前，需确认依赖和参数可复现。本轮没有运行任何新算法。

新验证集必须不同于已查看的 V1–V5、F1–F2。建议包含：纯相位误差、横向频率变化、鬼波、Q 正负误差、组合误差；合成数据至少保留独立物理生成，并增加二维传播场景。按地质模型/随机种子作为独立重复单位，不能将同一道集多个 Q 扰动当作独立样本。

建议待冻结门槛（尚未登记、不是实验结果）：复杂模型波形误差较 QScanRefit 降低至少 5%；Q 正确场景平均退化不超过 2%；弱反射、断层、尖灭各项平均退化不超过 2%；Q 扰动最坏误差不高于 QScanRefit。报告各模型分布、尾部失败率、相对振幅误差、频谱形状及计算耗时，不能只报告总体均值。全部参数与门槛在生成最终验证资料前冻结。

## 检索记录与边界

日期：2026-10-24；范围重点 2025-01-01 至检索日。网页搜索及出版社/作者资料核验；主要查询：
- seismic attenuation compensation 2026 Q wavelet joint estimation
- seismic blind nonstationary deconvolution wavelet attenuation 2025 2026
- seismic Q compensation physics informed 2025 unsupervised
- Nonstationary blind deconvolution using spectral constraints
- seismic Q 2026 joint wavelet estimation
- gxag053 DBSCAN

按 DOI 与标题去重。纳入上述 6 项，理由为可对应当前瓶颈、比较方法或工程实现；排除地震事件判别、地壳区域衰减、无直接迁移价值的一般去噪。另检出全波场反射/Q 联合估计会议工作，但对现有后处理代码改动过大，未列入近期方案。搜索引擎不提供可复现的完整命中总数，故不报告 PRISMA 数量；未执行 S2 去重，以 DOI/标题人工核对替代。来源以期刊为主，尚无跨数据集独立复现证据。

需纠正旧报告中可能过强的推断：一次 BAE 实现失败不等于整个 BAE 方法家族无效；QScan 的有效因子不等于真实物理 Q；联合子波与 Q 估计也不会自动消除可识别性问题。
