# 输出波形域联合事件收缩：代数与边界

本目录是固定 Q 与既有支持的机制 pilot。主候选与唯一 eta 网格见 PLAN_ZH.md。联合 Gram、相关噪声或凸二次收缩已有一般先例；查新范围与来源见相邻项目 innovation_direction_20261025/LITERATURE_REVIEW_ZH.md（包括 Pesquet 等 2009 §VI-B 式43）。本轮的可检验对象仅是本事件分组、实际 ridge 条件协方差和既有相干调制的具体组合。

## 上游实际坐标

上游 round3_20261024/run_combo.m 内嵌 refit_state，并非独立 refit_state.m。它调用 reflectivity_refit_prepared.m，再调用 reflectivity_soft_prepared.m；qscan_predictor 的 full_refit.m 也调用同一固定支持重拟合实现。原代码快照与提交/哈希均保存在 inputs/source_snapshots 和 inputs/source_manifest.json。

令 B 为衰减观测字典，每列归一化，scale 为归一化前列范数；W 是源波形字典。给定本道支持 S，C=B[:,S]，C=U diag(z) V'。上游定义

ρ=max(10^-8 z_max², (z_max²−10^4 z_min²)/(10^4−1))；
K=V diag(z/(z²+ρ)) U'；a=K d；r_S=a/scale_S。

因此不能用 OLS 的 (C'C)^(-1) 代替。若支持外生固定、观测噪声具有协方差 σ²I，则 a 的条件传播噪声协方差为 σ²H，H=KK'=V diag(z²/(z²+ρ)²)V'。本实现直接重建这个 ridge H，并用已归档反射率 a 检查正规方程 C'(Ca−d)+ρa=0，未重新发现支持。

设 P=W[:,S]/scale_S，组选择对角矩阵 J_k，Y_k=P J_k a。组 k 的噪声传播算子 A_k=P J_k K，所以

G_kl=Y_k'Y_l，N_kl=σ² tr(A_k A_l')=σ² tr(P_k H_kl P_l')。

代码以 H 与 P'P 的 Schur 乘积和分组指示矩阵 Z 计算 N=σ² Z'[(P'P)∘H]Z，与上述 trace 完全一致。上游 V_k 正是同一式的对角，E_k=Y_k'Y_k；数值检查必须先满足这两个等式。代码允许的差异仅为浮点舍入，并以对角合同缩放把 G/N 对角准确对齐原 E/V，保留 PSD；会报告对齐前相对误差，绝不把不同噪声尺度混在一起。

## 凸代理与对角退化

G 是向量 Y_k 的 Gram，N 是算子 σA_k 在 Frobenius 内积下的 Gram，所以两者半正定。D=diag(sqrt(λ_k))，λ_k=.5+15.5(1−c_k)²；原 fallback 组 λ=1。对角/全矩阵插值 M_eta=(1−eta)diag(M)+eta M 保持半正定。

盒约束目标为 F(g)=(g−1)'G_eta(g−1)+g'D N_eta Dg，0≤g≤1。H_q=G_eta+D N_eta D 半正定，目标连续且可行集紧，故最小值存在，解不保证唯一。梯度为 2[H_q g−G_eta 1]。盒 KKT：内部梯度为0，g=0处梯度≥0，g=1处梯度≤0。代码在缩放后的半梯度上报告归一化投影 KKT 残差，并在 L-BFGS-B 后做确定性坐标精修。

eta=0 时，目标逐组分离为 E_k(g_k−1)²+λ_k V_k g_k²，精确解 g_k=E_k/(E_k+λ_k V_k)；分母为0时沿用原 gain=1。这严格退化为同一 Python 实现的 Standard，含空支持与 fallback。eta=1 第一项是完整合成波形改变量 ||Σ_k(g_k−1)Y_k||²；0<eta<1 混合了此项与组独立能量代理。仅 G 耦合保留原 V 对角；仅 N 耦合保留原 E 对角，均使用主候选选定的同一 eta，不另选参。

## 为什么不是无偏风险

即使支持外生、g固定、D=I，记无噪声事件的 Gram 为 G_s，也有 E[G]=G_s+N。固定 g 的真实输出 MSE 为 (g−1)'G_s(g−1)+g'Ng。所用完整代理的期望则额外包含 (g−1)'N(g−1)；因此不是无偏 MSE 或 SURE。一般 DND 也不是未经调制的实际输出噪声能量。实际 g、c 与支持依赖数据，进一步超出了上述固定量计算。

联合最优仅保证降低相应代理，不保证真实波形误差降低、不保证收缩后波形能量降低、不保证弱事件改善，也不保证非扩张。0≤g≤1 不代表原子相消关系得到保留。每道 g 只重组本道事件波形；相干仍通过邻道影响其增益。条件白噪声传播未涵盖支持选择、校准 Q 的模型偏差、子波偏差或跨道噪声相关。

正交组的 G/N 均为对角时，eta=1 与 eta=0 理应一致；数值检查包含此负对照。整体与区域真值只用于预登记开发选 eta 及冻结后的诊断，不能由代理下降直接推出性能机制。
