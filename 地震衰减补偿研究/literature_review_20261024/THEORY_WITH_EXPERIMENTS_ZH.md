# 结合 GitHub 实验的补偿可恢复性推导与符号说明

日期：2026-10-24。通过 GitHub 连接器直接读取 q1754399106-sudo/work 默认分支的 9 个报告、代码和数据文件。未重跑 MATLAB；重新核算了输出误差分解 CSV 的 8 行恒等式。本文是理论分析，不是新方法验证结果。

## 1. 实验依据与对前一版推导的修正

1. q_identifiable_cluster_20261023/RESULTS_ZH.md：QAwareCluster 整体误差较 FullCalRefit 退化 2.40%，较 QAwareRidge 退化 3.04%。可辨识组合模式与降秩已经试过，不是尚未探索的创新。
2. q_output_trust_20261023/RESULTS_ZH.md：输出信赖域整体改善 0.18%，较 QAwareRidge 退化约 0.04%，未通过。限制偏离参考解不等于限制偏离真值。
3. q_output_trust_20261023/analysis/previous_round_output_diagnosis.csv：8/8 行 MSEChange>0；CrossTerm+ChangeNormTerm=MSEChange 的最大绝对差 4.17e-17。这里只核算已保存指标，不是重新生成波形。
4. bae_q_20261024/analysis/gates.json：QScan 独立组最坏误差比值均值 0.18815，Q0 误差比值均值 0.98383，但最差 Q0 比值 1.16912。均值改善不保证每个场景改善。独立弱反射值为 null，不能把 H4=true 解释成有该项数据证据。
5. complex_fieldlike_20261023/RESULTS_ZH.md：纯净目标保留真实源相位、接收增益及静校正；Moderate 的 FullCalRefit 振幅投影增益 0.61674，相关系数 0.85802，说明较高相关性不等于无振幅偏差。投影增益不是每个事件统一衰减比例。
6. q_error_stability_20261023/RESULTS_ZH.md：弱反射与整体误差可反向变化；复杂模型 Q0 是标定基准，不是真实横向 Q 全部正确。
7. METHOD_ZH.md 已使用剖除振幅/相位模板的投影和校准敏感性矩阵；新方案不能将“投影消去干扰参数”本身重新声称原创。

## 2. 统一符号，避免与代码混淆

所有向量默认列向量；以下使用实数时域表示，频域复数形式需将转置换为共轭转置。

| 符号 | 含义 | 维度或单位 |
|---|---|---|
| n | 处理块中观测样本总数，单道为时间采样点数，多道可展平 | 正整数 |
| p | 反射系数未知数数量 | 正整数 |
| k | 干扰参数总数 | 正整数 |
| m | 目标输出样本总数 | 正整数 |
| y | 实际观测地震数据 | R^n |
| r | 反射系数向量 | R^p，物理反射系数通常无量纲 |
| Q_l | 第 l 层品质因子，越小衰减通常越强 | 正数、无量纲 |
| ell_l=log Q_l | 第 l 层 Q 的对数参数，与仓库局部导数定义一致 | 无量纲 |
| beta | 子波频率、相位等参数；必要时扩展到增益、静校正 | 各分量分别为 Hz、rad、无量纲、s 等 |
| theta=(ell,beta) | 所有干扰参数组成的向量 | R^k |
| A(theta) | 从反射系数到含衰减观测的正演矩阵 | R^(n×p) |
| e | 观测噪声及未显式建模的误差 | R^n |
| x | 希望恢复的无吸收目标波形 | R^m |
| B(theta) | 从反射系数到目标波形的映射 | R^(m×p) |
| B0 | 预先固定的参考子波映射 | R^(m×p) |
| Cn | 工作噪声协方差，正定时可定义逆平方根 | R^(n×n) |
| L=Cn^(-1/2) | 噪声白化算子；这里用 L 避免与代码的子波 W 冲突 | R^(n×n) |
| z=(r,theta) | 全部未知量 | R^(p+k) |
| z0 | 局部线性化的中心解 | R^(p+k) |
| h=delta z | 相对中心解的小扰动 | R^(p+k) |
| F(z)=L A(theta)r | 白化后的预测观测 | R^n |
| T(z)=B(theta)r | 目标输出函数；T 不是传播时间 | R^m |
| J=DF(z0) | 白化预测对全部未知量的一阶导数 | R^(n×(p+k)) |
| G=DT(z0) | 目标输出对全部未知量的一阶导数 | R^(m×(p+k)) |
| Jr,Jtheta | J 对反射系数、干扰参数的分块 | n×p、n×k |
| Gtheta | 输出对干扰参数的导数 | m×k |
| epsilon | 白化总误差的界，包含声明的误差来源 | 白化数据范数单位 |
| epsilon_loc | 局部扰动允许的观测变化预算 | 白化数据范数单位 |
| u | 从输出提取一个成分的固定权重向量 | R^m |
| b=G^T u | 该输出成分对全部参数的敏感性 | R^(p+k) |
| I | 对应空间的单位矩阵 | 依上下文 |
| dagger | Moore–Penrose 伪逆 | 运算符 |
| ker | 零空间，即被矩阵映射到零的向量集合 | 集合 |
| range | 列空间 | 集合 |
| ||v||_2 | 向量 v 的欧氏长度 | 范数 |
| v^T w | 向量内积；上标 T 表示转置 | 标量 |
| hat | 估计值，如 xhat | 标记 |
| star | 真实值或明确指定的最优值，正文分别注明 | 标记 |

代码对应：METHOD_ZH.md 中 d 是本文 y，W 是子波字典，B 是列归一化的观测字典；代码 x 是归一化反射系数，s 才是输出波形。本文 x 统一表示输出，不能与代码 x 直接替换。qtrust_refit.m 中 C=d.B(:,support) 是当前支持上的观测字典，H=d.W(:,support)./d.scale(support)' 是输出映射；它的 H 不是本文的信息矩阵。

## 3. 两种目标必须分开

观测模型：y=A(theta)r+e。

标准化输出：x=B0 r。所有候选使用相同参考子波，目标导数 G=[B0,0]。

保留真实子波的去吸收输出：x=B(theta)r。为贴近复杂模型评分，B 应保留目标所要求的源相位、接收增益与静校正，仅去掉吸收及约定要去除的干扰。此时 G=[B(theta0),Gtheta]，Gtheta_i=(partial B/partial theta_i)r0。若 B 与 Q 无关，其 Q 导数为零；子波等参数的导数一般不为零。

改变评分目标需要从生成要素重新构造目标，并公平重评全部方法。不能将已保存的旧目标指标直接当作新目标效果。

## 4. 输出可识别的局部条件

F(z0+h)≈F(z0)+Jh，T(z0+h)≈T(z0)+Gh。

在线性化、无额外活动约束的模型中，输出可识别当且仅当 ker(J) 包含于 ker(G)。

证明：两个参数扰动产生相同观测意味着 J(h1-h2)=0；它们输出相同要求 G(h1-h2)=0。对所有这样的差向量成立，恰好就是零空间包含关系。

该结论只保证局部线性模型的唯一性，不保证稳定性。即使没有精确零空间，接近零的观测响应仍可对应很大的输出变化。非线性、多模态和不同支持解需要另行检查。

固定 L1 支持只能给出支持条件下的结论；仓库旧方法明确不能恢复支持阶段漏掉的事件。应将残差候选补入或多支持分析作为覆盖诊断，不能把固定支持内的不确定性当成全部不确定性。

## 5. 纳入独立校准，避免把已有信息丢掉

若独立校准提供白化线性观测变化 K delta theta，则联合观测导数为

Jall = [ Jr, Jtheta ; 0, K ]。

K 的行数等于校准约束数量，列数为 k。若校准还有未知振幅等参数，应先剖除这些自由度或将其一并纳入完整未知向量。K 必须对应可靠校准模型；一个带先验和特征截断的协方差代理不能自动充当已校准的观测信息。

定义完整局部信息矩阵 H=Jall^T Jall，大小 (p+k)×(p+k)。

在标准化目标下剖除 theta，得到

Heff = Jr^T Jr - Jr^T Jtheta (Jtheta^T Jtheta + K^T K)^dagger Jtheta^T Jr。

Heff 大小 p×p，描述干扰参数可调整但受校准约束时，数据对反射系数的剩余约束。这个表达式来自最小化 ||Jr delta r+Jtheta delta theta||²+||K delta theta||²。没有校准时 K=0，退化为上一版的正交投影公式。

对于保留未知子波的输出，不能只用 Heff 和固定 B；因为输出自身随 theta 改变，须使用完整的 Jall 与 G。

## 6. 输出成分的局部不确定性

s=u^T x 表示一个指定输出成分，u 可以是单点采样权重或一个固定实值分析模板。其扰动 delta s=b^T h，其中 b=G^T u。

对局部扰动集合 h^T H h≤epsilon_loc²，定义

U_u = sup |b^T h|。

若 b 与 ker(H) 正交，则 U_u=epsilon_loc sqrt(b^T H^dagger b)；否则在没有额外参数界的线性模型里 U_u=∞。

证明：在 H 的正特征值子空间写 H=V diag(lambda_i)V^T，用加权 Cauchy–Schwarz 不等式得到 |b^T h|²≤(sum (v_i^T b)²/lambda_i)(h^T H h)，取对应比例方向达到等号。若 b 在零空间有非零投影，沿该方向放大 h 不增加观测变化却无界改变输出。

V 的列 v_i 是单位特征向量，lambda_i 是非负特征值；U_u 是该输出成分的允许变化上界。这里的无穷是无限制线性模型的结果，不代表物理参数可以无限变化；实际有界域需要约束求解。

若中心 z0 不是精确拟合，任一可接受候选与中心的预测差满足 ||F(z)-F(z0)||≤epsilon+||Ly-F(z0)||。线性化还需计入余项界。因此 epsilon_loc 不是自动等于 epsilon。候选间比较的观测差上界为 2epsilon。

## 7. 旧信赖域为什么不能保证改善

设 xref 是完整参考补偿，xtrue 是真值，delta=xnew-xref 是修改量，eref=xref-xtrue 是参考误差。则

||xnew-xtrue||²-||xref-xtrue||² = 2 eref^T delta + ||delta||²。

改善必须满足 eref^T delta < -||delta||²/2。

旧信赖域约束了输出变化 delta 的部分投影，却不知道真实误差 eref 的方向。小步长最多限制退化幅度，不能保证改变方向正确。CSV 的 8 行均为正误差变化与这个结论一致，但不证明每次降秩都会失败。

对应代码实际使用固定先导的 Q 误差白化、簇内 SVD、带轻岭的联合重拟合，以及参考输出投影超过 2 个近似标准差后的累计变化约束。它不是对所有物理可接受模型的最大输出差求解。

## 8. 面向输出的全局集合与最坏误差

定义 C(y)={z∈Z:||Ly-F(z)||≤epsilon，且符合独立校准约束}。Z 是公开声明的物理和结构约束域。X(y)={T(z):z∈C(y)} 是允许的输出集合。

D(y)=sup_{x1,x2∈X(y)}||x1-x2|| 是直径。任何估计器的最坏误差至少为 D(y)/2；最小包围球中心在非空紧致集合下存在，其半径 R 满足 D/2≤R≤D。

所有保证都要求真实模型/输出落在所建集合内。模型族失配、错误噪声界、遗漏支持都会破坏这个前提。有限扫描得到的集合通常偏小：已发现两个相距远的可行解可证明歧义存在；未找到相距远的解不能证明没有歧义。

## 9. Q 不确定性对应的简化频带边界

局部单分量振幅模型 Y(f)=exp(-pi f tau/Q) X(f)。Y 是观测频谱，X 是目标频谱，f 是频率(Hz)，tau 是该模型的传播时间(s)，pi 是圆周率。这里忽略噪声、相位和干涉，并假设源子波影响已处理。

给定 Q∈[Qmin,Qmax]，目标振幅 a=|X| 位于 [aL,aU]，其中 aL=|Y|exp(pi f tau/Qmax)，aU=|Y|exp(pi f tau/Qmin)。

最坏相对振幅误差最小的估计为 ahat=2aL aU/(aL+aU)，误差下界 eta_min=(aU-aL)/(aU+aL)=tanh[pi f tau(1/Qmin-1/Qmax)/2]。

若希望误差不超过 eta∈(0,1)，必须 f≤2 artanh(eta)/[pi tau(1/Qmin-1/Qmax)]。artanh 是反双曲正切，eta 为容许相对误差；Q 区间退化时该式只表明没有来自 Q 不确定性的频率限制，噪声等限制仍存在。

分层模型可改用吸收积分 t*=sum_l tau_l/Q_l；tau_l 为各层传播时间，t* 的单位为秒。如果 t*∈[t*min,t*max]，则 eta_min=tanh[pi f(t*max-t*min)/2]。这样更贴合仓库累计 dt/Q 的正演。

## 10. 前一版其他符号对照

q=1/Q：前一版使用的逆品质因子；本版改用 ell=logQ 以对应仓库代码。两者导数转换为 partial F/partial ell_l=-(1/Q_l) partial F/partial q_l，不可混用协方差尺度。

P_perp=I-Jtheta Jtheta^dagger：无校准、自由干扰参数情况下的正交补投影。H_eff：剖除干扰后的信息矩阵。D：输出集合直径。R_mm：最小最大误差半径。下标 mm 表示 minimax。

上一版收缩推导中的 M：正定参数椭球尺度矩阵，约束 h^T M^-1 h≤1，不默认表示真实后验协方差；sigma：白噪声标准差；a：线性估计器权重，不是本版频谱振幅；rho_i：第 i 个参数模态允许幅度；s_i：该模态观测灵敏度；g_i：相对于直接逆解的收缩系数，范围 0 到 1。lambda_i：信息矩阵特征值。delta：扰动或相对误差时必须按所在公式说明，不能跨式默认同义。

## 11. 可发展的贡献与下一步

已有实验已经覆盖 Q 白化、簇模式、输出信赖域；新增贡献不能只是这些组件的重命名。值得检验的是：用完整观测与校准联合约束、正确的目标导数，预测“哪些输出不能可靠恢复”；再将这种边界用于指定保幅误差下的补偿。

第一轮只做诊断：固定旧算子和数据，分别计算仅 Q、Q+主频/相位、加入独立校准后的敏感性；检查预测不稳定的输出成分是否对应已有失真。不要立刻根据小奇异值删模式。线性化可用有限差分对照，真实误差只用于事后评分。上述新诊断尚未运行。

经典基础：Donoho (1994), Statistical Estimation and Optimal Recovery, https://doi.org/10.1214/aos/1176325367 。零空间判据、投影、Schur 补和最小最大恢复均有已有理论基础，本笔记不宣称这些基本工具原创。

## GitHub 来源

https://github.com/q1754399106-sudo/work/blob/main/q_identifiable_cluster_20261023/RESULTS_ZH.md
https://github.com/q1754399106-sudo/work/blob/main/q_identifiable_cluster_20261023/METHOD_ZH.md
https://github.com/q1754399106-sudo/work/blob/main/q_output_trust_20261023/RESULTS_ZH.md
https://github.com/q1754399106-sudo/work/blob/main/q_output_trust_20261023/core/qtrust_refit.m
https://github.com/q1754399106-sudo/work/blob/main/q_output_trust_20261023/analysis/previous_round_output_diagnosis.csv
https://github.com/q1754399106-sudo/work/blob/main/bae_q_20261024/analysis/gates.json
https://github.com/q1754399106-sudo/work/blob/main/q_error_stability_20261023/RESULTS_ZH.md
https://github.com/q1754399106-sudo/work/blob/main/complex_fieldlike_20261023/RESULTS_ZH.md
https://github.com/q1754399106-sudo/work/blob/main/JOINT_REGULARIZATION_DIRECTIONS_20261022.md
