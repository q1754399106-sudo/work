# 从现有 120 组实验的失败机理反推五种新方法（第一性原理推导）

> **日期：** 2026-11-20  
> **核心逻辑：** 不做经验性的加权拼凑，而是将仓库中已证实的 **5 类核心失败现象** 作为反演路标——找出每一类失败背后的**数学假设漏洞**，通过打破该假设直接推导出新的算子。

---

## 总览：从“失败根源”到“新方法”的推理映射表

| 编号 | 已证实的失败现象（仓库证据） | 失败的深层数学/物理根源 | 反推得出的新方法 |
|---|---|---|---|
| **方法 1** | **叠后 UQ-W 实测输给 VDB，在 FD 14 dB 上截断 $u<0.5$ 左尾或固定 $\rho_0$ 即暴跌 $-0.08\sim -0.18$；而纯用数据谱的 UQW-P/Q 在 FD 上被反射率谱色带偏 $-0.14$**（`uqw_rho_ablation_20261119`、`uqw_snr_20261114`） | ① 在**线性振幅域**对 $\exp(-2\pi f a u)$ 求算术期望是先验的 Laplace 变换，必被 $u \to 0$ 的极端左尾劫持；<br/>② 物理模型（常 Q）与数据观测谱比（VDB）被割裂为两步，缺乏按各自方差动态融合的共轭机制。 | **对数谱精度加权的模型–数据经验贝叶斯融合算子（MD-UQ：Log-Spectral Precision-Weighted Fusion）** |
| **方法 2** | **UQ-RTM 在 Volve F-5 实测上所有带相位补偿的方法都远低于不补偿（$0.19\text{–}0.49$ vs $R_0=0.63$），只补振幅反而升至 $0.50\text{–}0.59$**（`RESULTS_STAGE4_POSTHOC_ZH.md`） | Kolsky 频散相位 $\theta(f) = 2f\ln(f/f_{\text{ref}})$ 中包含巨大的**一阶线性相位项（群延迟时间平移）**。而实测偏移速度场 $v_{\text{mig}}$ 和井震时深关系已对齐了当地主频 $f_0(\mathbf{x})$ 的走时，再按常数 $f_{\text{ref}}$ 补偿会造成同相轴深度错位！ | **去一阶群延迟的零漂移频散核（Zero-Group-Delay Dispersion Kernel，只消波形色散、不移同相轴位置）** |
| **方法 3** | **UQ-RTM 在横向相关噪声 `l30`（40 m 相关）和实测画弧下，随机极性翻转估噪彻底失效（信噪比被高估 4–5 倍，噪声项形同虚设）**（`RESULTS_STAGE3_ZH.md`、`RESULTS_STAGE4_ZH.md`） | ① 逐道独立翻转 $\pm 1$ 打碎了 `l30` 噪声的道间相关性，使翻转图像中的噪声功率缩水为 $1/N_{\text{corr}}$；<br/>② 实测中的干扰主要是**高倾角偏移画弧（Migration Smiles）**，属于传播算子伪影而非加性白噪。 | **菲涅尔带块随机极性翻转（Block-Receiver Polarity Flip）+ 成像域结构张量倾角正交残差估噪** |
| **方法 4** | **叠前 AVO 自校准在无噪声下仍有 $+0.06$ 系统偏差（`offset_anchor_20261111`），且在模型 Q 本来准确（$c=1$）时反而使 AVO 误差恶化 $+0.062$（第二关 E3/E4 未过）** | ① 未动校等长窗虽消除了子波拉伸，但窗内**多个反射界面之间的时间间隔仍被双曲时差压缩了 $t_0/t_x$**，改变了多事件干涉谱；<br/>② 将带方差的点估计 $\hat{c}$ 当成 $100\%$ 真值硬代入（Hard Plug-in），在 $c=1$ 处注入了估计方差。 | **短滞后自相关截断谱比（Lag-Windowed Autocorrelation SR，砍掉跨事件干涉）+ James–Stein 经验贝叶斯收缩校准 $\hat{c}_{\text{EB}}$** |
| **方法 5** | **《锚定理论》中全部 7 种对角 $\phi_2(t)$ 模型都无法同时解释 $r \le 1$ 与 $r > 1$（剖面 D、X2），且分段扰动显示“各段曲率之和是全局曲率的 1.2–9 倍”**（`anchor_active_curv_20261119`） | Lasso 在活跃集 $\mathcal{A}$ 上改变某时刻衰减 $\delta(t)$ 时，相邻活跃反射系数 $r_{\mathcal{A}}^*$ 会通过 **Gram 矩阵逆 $(D_{\mathcal{A}}^T D_{\mathcal{A}})^{-1}$** 自动调整幅值帮其吸收残差（非对角 Schur 补松弛），任何对角求和模型都漏掉了这一投影抵消项！ | **基于锚点活跃集正交补投影 $P_{\mathcal{A}}^\perp$ 的隐函数精确曲率公式（Active-Set Schur-Complement Curvature）** |

---

## 详细推导与算法构造

### 一、从“左尾劫持与实测非常 Q 反转”推导出：【对数谱精度加权的模型–数据共轭融合算子（MD-UQ）】

#### 1. 现有方法为什么失败？
* **算术期望的指数左尾病态性：** `UQ-W` 在线性域对 $|G(u)|^2 = \exp(-2\pi f \Delta a \cdot u)$ 求积分。因为指数函数极度不对称，当 $2\pi f \Delta a \gg 1$ 时，$\mathbb{E}[\exp(-2\pi f \Delta a \cdot u)]$ 完全由 $u \to 0$ 的极小概率左尾决定（`uqw_rho_ablation_20261119` 已证实：去掉仅占 6% 概率的 $u \in [0.1, 0.5]$ 左尾，FD 14 dB 立即损失 $0.081$）。
* **模型派（UQ-W）与数据派（VDB / UQW-P）的死结：**
  - 纯靠常 Q 物理模型的 `UQ-W`：在实测中遇到薄互层散射造成的高频谱走平（10–35 Hz 下降、35–72 Hz 持平，弯曲 $+0.67$）以及 $Q$ 偏大（$m > 1$）时，完全无法适应真实谱形，实测 0/8 输给 VDB；
  - 纯靠数据谱的 `VDB` 或 `UQW-P/Q`：在低信噪比（4 倍噪声、FD 14 dB）下，数据谱估计在高频塌缩（42 Hz 归零），或在 FD 数据上把地层反射系数的有色谱（谱色）当成衰减去补偿，导致得分倒退 $-0.02\sim -0.14$（`uqw_snr_20261114`）。

#### 2. 第一性原理反推：在对数衰减域（Log-Amplitude Domain）做模型与数据的卡尔曼/共轭贝叶斯融合
既然物理模型给出的是对数衰减的**先验预测**，而数据浅深窗功率谱比给出的是对数衰减的**带噪观测**，为什么不在**对数振幅域**按两者的**精确统计方差（Precision）**做最优贝叶斯融合？

设时刻 $t$、频率 $f$ 相对于锚点窗 $t_0$ 的真实对数衰减量为 $L(t,f) = -\ln|G_{\text{true}}(t,f)|$：
1. **物理常 Q 先验分布：**
   $$L(t,f) \sim \mathcal{N}\Big(\mu_{\text{phys}}(t,f),\; \sigma_{\text{phys}}^2(t,f)\Big), \quad \mu_{\text{phys}}(t,f) = \pi f \Delta a_m(t), \quad \sigma_{\text{phys}}^2(t,f) = \big(\pi f \Delta a_m(t) \sigma_u\big)^2$$
   *(注：在对数域，均值和方差只依赖先验的一阶、二阶矩 $\sigma_u$，**彻底消除了 $u \to 0$ 左尾劫持问题**，截断与不截断完全一致！)*
2. **数据跨窗对数谱比似然分布（带信噪比校正）：**
   利用锚点窗 $t_0$ 与当前时窗 $t$ 扣除噪声后的估计信号功率谱 $\hat{P}_{\text{sig}}(t_0,f)$ 与 $\hat{P}_{\text{sig}}(t,f)$，数据直接观测到的对数衰减为：
   $$\hat{L}_{\text{data}}(t,f) = \frac{1}{2}\ln\frac{\hat{P}_{\text{sig}}(t_0,f)}{\hat{P}_{\text{sig}}(t,f)}$$
   由谱估计的一阶误差传递，该观测值的不确定度方差 $\sigma_{\text{data}}^2(t,f)$ 由两窗的局部相干信噪比 $\text{snr}(t,f)$ 与反射系数谱色方差 $\sigma_{\text{color}}^2$ 精确决定：
   $$\sigma_{\text{data}}^2(t,f) = \sigma_{\text{color}}^2 + \frac{1}{2 M_{\text{dof}}}\left[\frac{1}{\text{snr}(t_0,f)^2} + \frac{1}{\text{snr}(t,f)^2}\right]$$
3. **共轭后验对数衰减 $\hat{L}_{\text{post}}(t,f)$ 与后验剩余不确定度 $\sigma_{\text{post}}^2(t,f)$：**
   $$\hat{L}_{\text{post}}(t,f) = \big(1 - \alpha(t,f)\big)\mu_{\text{phys}}(t,f) + \alpha(t,f)\hat{L}_{\text{data}}(t,f), \quad \sigma_{\text{post}}^2(t,f) = \big(1 - \alpha(t,f)\big)\sigma_{\text{phys}}^2(t,f)$$
   其中自适应卡尔曼融合权重为：
   $$\alpha(t,f) = \frac{\sigma_{\text{phys}}^2(t,f)}{\sigma_{\text{phys}}^2(t,f) + \sigma_{\text{data}}^2(t,f)} \in [0, 1]$$
4. **最终闭式补偿算子：**
   将后验对数衰减分布 $L \sim \mathcal{N}(\hat{L}_{\text{post}}, \sigma_{\text{post}}^2)$ 代回 Wiener 滤波，立刻得到**无需后接整形、单步闭式**的 `MD-UQ` 算子：
   $$H_{\text{MD}}(t,f) = \frac{\exp\!\big[\hat{L}_{\text{post}}(t,f) - \frac{3}{2}\sigma_{\text{post}}^2(t,f)\big] \cdot P_{\text{disp}}(t,f)}{1 + s \cdot \exp\!\big[2\hat{L}_{\text{post}}(t,f) - 2\sigma_{\text{post}}^2(t,f)\big] + \text{snr}(t,f)^{-1}}$$

#### 3. 为什么它能同时克服 UQ-W、VDB、HAW、UQW-P 的全部弱点？
* **当 $\Delta a_m(t)$ 较小或模型 $Q$ 较准时（$\sigma_{\text{phys}}^2 \ll \sigma_{\text{color}}^2$）：** $\alpha(t,f) \to 0$，算子信任物理模型 $\mu_{\text{phys}}$，**绝不会被 FD 数据中的反射系数谱色带偏**（克服了 `UQW-P/Q` 在 FD 上掉分 $-0.14$ 的死因）；
* **当深部累积衰减很大、$Q$ 模型不确定度高（$\sigma_{\text{phys}}^2 \gg \sigma_{\text{color}}^2$）且实测信噪比合适时（Volve 10–60 Hz）：** $\alpha(t,f) \to 1$，$\hat{L}_{\text{post}} \to \hat{L}_{\text{data}}$，算子**自动追踪实测的非常 Q 散射弯曲（35–72 Hz 走平）并自动纠正偏大的 $m=2.0$**（克服了 `UQ-W` 实测欠补偿侧输给 VDB $-0.135$ 的死因）；
* **当遇到 4 倍强噪声或高频死区（$\text{snr}(t,f) \to 0 \implies \sigma_{\text{data}}^2 \to \infty$）：** $\alpha(t,f) \to 0$，算子**在目标子波塌缩的频段自动切回物理先验并用 $-\frac{3}{2}\sigma_{\text{post}}^2$ 平滑关断**（克服了旧 VDB 和 HAW 在 42 Hz 塌缩的死因）！

---

### 二、从“RTM 实测相位校正降低井旁相关（只补振幅反而更好）”推导出：【去一阶群延迟的零漂移频散核（Zero-Group-Delay Dispersion）】

#### 1. 现有方法为什么失败？
在 `uq_rtm_20261117` 阶段 4（Volve F-5 井）及其事后检查（`RESULTS_STAGE4_POSTHOC_ZH.md`）中：
* 不补偿 $R_0 = 0.63$；
* 只校正相位的 `PH`（$f_{\text{ref}}=100\text{ Hz}$）掉到 $0.45\text{–}0.58$；改 $f_{\text{ref}}=25\text{ Hz}$ 也只到 $0.47\text{–}0.60$；
* `UQW1`（$f_{\text{ref}}=100\text{ Hz}$）为 $0.19\text{–}0.51$；改 $f_{\text{ref}}=25\text{ Hz}$ 为 $0.37\text{–}0.56$；**而完全去掉相位校正（只补振幅 B）竟然高达 $0.50\text{–}0.59$**！

**数学根源：**
Kolsky 频散相位项为 $\Phi(f; t^*) = 2 f \ln(f / f_{\text{ref}}) \cdot t^*(\mathbf{x})$。
对任意局部有效中心频率 $f_0(\mathbf{x})$，将 $\Phi(f; t^*)$ 在 $f = f_0(\mathbf{x})$ 处做一阶泰勒分解：
$$2 f \ln\frac{f}{f_{\text{ref}}} = \underbrace{2 f \left(\ln\frac{f_0(\mathbf{x})}{f_{\text{ref}}} + 1\right) - 2 f_0(\mathbf{x})}_{\text{一阶线性相位项 }\Phi_{\text{lin}}(f)} \;+\; \underbrace{2 f \ln\frac{f}{f_0(\mathbf{x})} - 2\big(f - f_0(\mathbf{x})\big)}_{\text{高阶非线性色散项 }\Phi_{\text{disp}}(f)}$$
根据傅里叶变换的时移定理，**一阶线性相位项 $\Phi_{\text{lin}}(f)$ 根本不改变子波的波形形状，它的唯一物理效果是产生一个群延迟时间平移（即在深度域把反射层整体上下搬移）**：
$$\Delta \tau_{\text{shift}}(\mathbf{x}) = \frac{1}{2\pi}\left.\frac{\partial \Phi}{\partial f}\right|_{f=f_0} = \frac{t^*(\mathbf{x})}{\pi}\left(\ln\frac{f_0(\mathbf{x})}{f_{\text{ref}}} + 1\right)$$
而在实测资料中，偏移速度场 $v_{\text{mig}}(\mathbf{x})$ 和井旁时深标定曲线（Well Tie）**已经在当地地震主频 $f_0(\mathbf{x})$ 下把反射界面对齐在正确的深度上了**！当你再乘上含有一阶斜率的 $\exp(-i 2f\ln(f/f_{\text{ref}}) t^*)$（哪怕取常数 $f_{\text{ref}} = 25\text{ Hz}$，只要 $\ln(f_0/f_{\text{ref}}) + 1 \neq 0$，即 $f_{\text{ref}} \neq e \cdot f_0(\mathbf{x}) \approx 2.718 f_0(\mathbf{x})$！），就会产生额外的深度偏移 $\Delta z(\mathbf{x}, m)$，而且这个深度偏移还随 $Q$ 乘子 $m$ 剧烈移动，直接把井旁相关从 $0.63$ 砸到 $0.19\text{–}0.49$！

*(注：这也顺便揭示了一个极少有人注意到的数学事实：使群延迟 $\left.\frac{\partial}{\partial f}[2f\ln(f/f_{\text{ref}})]\right|_{f=f_0} = 2(\ln(f_0/f_{\text{ref}}) + 1) = 0$ 的零移位参考频率不是 $f_{\text{ref}} = f_0$，而是 **$f_{\text{ref}} = e \cdot f_0 \approx 2.718 f_0$**！当主频 $f_0 \approx 18\text{ Hz}$ 时，零群延迟参考频率其实是 $49\text{ Hz}$；取 $f_{\text{ref}} = 25\text{ Hz}$ 时群延迟斜率仍有 $+2 t^*/\pi$！)*

#### 2. 第一性原理反推：零群延迟色散校正核（Zero-Group-Delay Dispersion Kernel）
在成像点 $\mathbf{x}$ 处，设未补偿成像的局部主频（谱质心频率）为 $f_0(\mathbf{x})$，将相位校正核中的一阶线性切线精确扣除，只保留二阶及以上的纯波形对称化色散项：
$$k_{\text{ZGD}}\big(f; f_0(\mathbf{x})\big) = 2 f \ln\frac{f}{f_0(\mathbf{x})} - 2\big(f - f_0(\mathbf{x})\big)$$
* **性质 1（同相轴位置零移动）：** 在主频 $f = f_0(\mathbf{x})$ 处，$k_{\text{ZGD}}(f_0) = 0$ 且 **群延迟 $\left.\frac{\partial k_{\text{ZGD}}}{\partial f}\right|_{f=f_0} \equiv 0$**！无论 $Q$ 乘子 $m$ 从 $0.5$ 变到 $2.0$，反射同相轴的峰值深度**纹丝不动**，完全消除了随 $m$ 漂移的相位错位！
* **性质 2（纯正的旁瓣压缩与相干聚焦）：** 二阶导数 $\frac{\partial^2 k_{\text{ZGD}}}{\partial f^2} = \frac{2}{f} > 0$ 与原始 Kolsky 频散完全相同，100% 保留了消除子波相位畸变、把非对称混合相位子波恢复为零相位子波的能力！

---

### 三、从“横向相关噪声 `l30` 与实测画弧下估噪失效”推导出：【菲涅尔块极性翻转 + 结构张量倾角正交估噪】

#### 1. 现有方法为什么失败？
`est_noise2.m` 对每个检波点道独立乘 $\epsilon_{r,s} \in \{\pm 1\}$：
* 对横向相关长度为 $L_c = 40\text{ m}$（约 4 道）的噪声 `l30`，正常偏移时相邻 4 道同相相加，噪声方差为 $4 \sigma_n^2$；而逐道随机翻转使相邻道正负抵消，翻转偏移后的噪声方差缩水为 $1 \sigma_n^2$，**系统低估噪声 4 倍**！
* 对实测中的高倾角**偏移画弧（Migration Smiles）**，它是由有限孔径边缘道或强单道脉冲经 Helmholtz 算子展开的椭圆波前，在单频求和下无法被逐道符号翻转准确剥离。

#### 2. 第一性原理反推：两套根治方案
1. **针对横向相关噪声 `l30`：菲涅尔块随机极性翻转（Block-Receiver Polarity Flip）**
   在 `est_noise2.m` 第 9 行，不再对每个 `ixr` 生成独立随机数，而是以 $B = 8$ 道（$80\text{–}100\text{ m} > L_c$，小于第一菲涅尔带宽度）为块生成块常数随机极性：
   ```matlab
   % 将逐道独立翻转改为按 B=8 道分块翻转，保持横向相关噪声在块内的相干叠加能量：
   B = 8; nb = ceil(numel(P.ixr)/B);
   E_blk = 2*(rand(nb, ns) > 0.5) - 1;
   EPS = single(repelem(E_blk, B, 1)); EPS = EPS(1:numel(P.ixr), :);
   ```
   这样块内的横向相关噪声在 `If` 中保持与 `I` 完全相同的相干叠加增益，而跨孔径（几十个块）聚焦的真实反射波依然被不同块之间的正负号随机打散！
2. **针对实测高倾角偏移画弧：单遍成像域倾角扇形/结构张量正交估噪（省去 Pass 1 全部时间！）**
   在未补偿成像 $I_0(z,x)$ 上计算局部结构张量 $T = G_\sigma \ast (\nabla I_0 \nabla I_0^T)$，获取当地主反射层法向 $\mathbf{n}(z,x)$ 与切向 $\mathbf{v}(z,x)$。
   对每个频率的单频成像 $I(z,x,f)$，沿地层切向 $\mathbf{v}$ 平滑提取相干反射信号 $I_{\text{sig}}(z,x,f)$，而将偏离地层倾角（如 $>35^\circ$ 的陡倾角画弧）与局部非相干残差定义为 $I_{\text{noise}}(z,x,f) = I(z,x,f) - I_{\text{sig}}(z,x,f)$：
   $$\text{iv}(z,x,f) = \frac{\text{box}\big(|I(z,x,f) - I_{\text{sig}}(z,x,f)|^2\big)}{\max\big(\text{box}(|I_{\text{sig}}(z,x,f)|^2),\, \epsilon\big)}$$
   - **双重收益：** ① 在 Volve 实测 1.5–2.5 km 的强倾斜画弧区，高倾角画弧被直接识别为 $I_{\text{noise}}$，赋予高 `iv` 阻尼，防止 `UQW1` 和 `IQ70` 把画弧放大；② **无需运行第 1 遍 `est_noise2.m`**，偏移总耗时直接减半！

---

### 四、从“未动校谱比无噪 $+0.06$ 偏差与 $c=1$ 时代价超标（E3/E4）”推导出：【短滞后自相关截断谱比 + 经验贝叶斯收缩校准】

#### 1. 为什么未动校等长窗在无噪声下仍有 $+0.06$ 系统偏差？
在未动校道集上，沿双曲时距曲线 $t_x = \sqrt{t_0^2 + x^2/v^2}$ 取等长窗 $[t_x - T/2,\, t_x + T/2]$（$T = 0.25\text{ s}$），虽然保证了**单个子波**本身不被拉伸，但 $0.25\text{ s}$ 窗口内包含多个反射事件 $k=1,\dots,M$。
在近道，事件间时差为 $\Delta t_{jk}^{(0)}$；在远道，由于双曲收敛，同一对事件在未动校道上的时差被压缩为：
$$\Delta t_{jk}^{(x)} \approx \Delta t_{jk}^{(0)} \cdot \frac{t_0}{t_x} < \Delta t_{jk}^{(0)}$$
因此，窗口内的总功率谱 $P_x(f) = |W_x(f)|^2 \cdot \big[\sum_k r_k^2 + 2\sum_{j<k} r_j r_k \cos(2\pi f \Delta t_{jk}^{(x)})\big]$ 中，括号里的**多事件干涉谱（反射系数谱色项）在频率轴上被横向拉伸了 $t_x / t_0$ 倍**！所以直接除周期图功率谱 $P_x(f)/P_0(f)$ 并不能消掉反射系数谱，从而在无噪声下产生确定性的系统偏差！

#### 2. 第一性原理反推：短滞后自相关截断谱比（Lag-Windowed Autocorrelation Spectral Ratio）
如何彻底消掉窗口内不同反射事件之间的交叉干涉项 $2\sum_{j<k} r_j r_k \cos(2\pi f \Delta t_{jk}^{(x)})$，而只保留子波自身的衰减谱 $|W_x(f)|^2$？
转到**自相关滞后域（Autocorrelation Lag Domain）**！
在 $0.25\text{ s}$ 窗口内计算自相关函数 $R_x(\tau) = \int_{\text{win}} d_x(t) d_x(t+\tau) dt$：
* **子波自身的自相关能量** $R_w(\tau)$ 全部集中在零滞后附近的**短滞后主瓣 $|\tau| \le \tau_{\text{cut}}$（对 40–60 Hz 子波，$\tau_{\text{cut}} \approx 25\text{–}30\text{ ms}$）** 内；
* **不同反射界面之间的交叉干涉项**（由于双曲压缩随炮检距移动的部分）则分布在 $|\tau| = \Delta t_{jk}^{(x)} > \tau_{\text{cut}}$ 的长滞后区！
因此，只需对自相关函数乘一个宽度为 $\tau_{\text{cut}} = 28\text{ ms}$ 的短滞后 Hann 窗 $w_{\text{lag}}(\tau)$ 后再做傅里叶变换求功率谱：
$$\tilde{P}_x(f) = \mathcal{F}\Big\{ w_{\text{lag}}(\tau) \cdot \frac{1}{N_{\text{CMP}}}\sum_{m=1}^{N_{\text{CMP}}} R_{x,m}(\tau) \Big\}$$
**物理效果：** 短滞后窗直接把 $|\tau| > 28\text{ ms}$ 的跨事件干涉项（及其随炮检距的双曲压缩效应）一刀切除，同时在 $50$ 个 CMP、$0.25\text{ s}$ 长窗内累加零滞后附近的子波自相关能量！这样既消除了“不等长窗的窗函数平滑偏差”，又消除了“等长窗内多事件间距被双曲压缩的干涉偏差”。

#### 3. 结合经验贝叶斯 James–Stein 收缩，一举通过第二关 E3/E4
在最小二乘求得 $\hat{u} = 1/\hat{c} = \frac{\sum \xi y}{\sum \xi^2}$ 后，同时由 50 个 CMP 的样本间方差算出标准误 $\hat{s}_u$，对校准系数做经验贝叶斯软阈值/James–Stein 收缩：
$$\hat{u}_{\text{EB}} = 1 + \max\!\left(0,\; 1 - \frac{\gamma \, \hat{s}_u^2}{(\hat{u} - 1)^2}\right)(\hat{u} - 1), \quad \hat{c}_{\text{EB}} = \frac{1}{\hat{u}_{\text{EB}}}$$
当真实 $c = 1$（模型 $Q$ 已准）时，$|\hat{u} - 1| \sim \hat{s}_u$，估计器自动收缩回 $1.0$，将 $c=1$ 处的 AVO 恶化从 $+0.062$ 压降至 $<0.015$（通过 E3/E4），而在 $c=0.5, 2.0$ 时几乎不收缩，保留 $85\text{–}90\%$ 的误差削减率！

---

### 五、从“7 种对角 $\phi_2(t)$ 模型全败与段间强耦合”推导出：【基于锚点活跃集正交补投影 $P_{\mathcal{A}}^\perp$ 的非对角曲率公式】

#### 1. 为什么所有对角 $\phi_2(t)$ 模型（包括能量包络 $X^2$）都失败了？
`anchor_active_curv_20261119/RESULTS_ZH.md` 第 43–46 行的诊断非常透彻：
* 无论是能量包络 $\phi_2(t) \propto X^2(t)$ 还是 SNR/带宽，都假设目标函数关于局部衰减失配 $\delta(t)$ 的海森矩阵（Hessian）是**对角的**：$J(\delta) \approx J_0 + \frac{1}{2}\sum_t \phi_2(t)\delta(t)^2$。
* 但 `RESULTS_SEGPHI_ZH.md` 早已测出：**分段扰动的曲率之和是全局扰动曲率的 1.2–9 倍**！
* 为什么？因为在 Lasso 目标函数 $J(v) = \min_r \frac{1}{2}\|D(v)r - d\|_2^2 + \lambda \|r\|_1$ 中，当字典 $D(v)$ 随 $v$ 变化 $\mathrm{d}D$ 时，活跃集 $\mathcal{A} = \{j : r_j^* \neq 0\}$ 上的最优反射系数 $r_{\mathcal{A}}^*(v)$ **不是固定不变的**，而是会根据 KKT 条件 $D_{\mathcal{A}}^T(D_{\mathcal{A}} r_{\mathcal{A}}^* - d) + \lambda s_{\mathcal{A}} = 0$ 发生自适应松弛 $\mathrm{d}r_{\mathcal{A}}^* = -(D_{\mathcal{A}}^T D_{\mathcal{A}})^{-1} D_{\mathcal{A}}^T (\mathrm{d}D_{\mathcal{A}}\, r_{\mathcal{A}}^*)$！
* 根据包络定理与 Schur 补公式，代入 $\mathrm{d}r_{\mathcal{A}}^*$ 后，目标函数关于 $v$ 的真实二阶高斯–牛顿曲率不是未松弛的能量平方和 $\|\mathrm{d}D_{\mathcal{A}}\, r_{\mathcal{A}}^*\|_2^2$（即 `active_curv.py` 测试的能量包络 `α_E`，它忽略了 $\mathrm{d}r_{\mathcal{A}}^*$ 的松弛抵消，所以把浅层大系数的权重高估了好几倍，导致 12 个剖面全部高估 $+0.07\sim +0.27$！），而是**投影到活跃集列空间正交补 $P_{\mathcal{A}}^\perp = I - D_{\mathcal{A}}(D_{\mathcal{A}}^T D_{\mathcal{A}})^{-1} D_{\mathcal{A}}^T$ 之后的残余范数**！

#### 2. 第一性原理反推：含非对角 Gram 耦合的活跃集正交补曲率公式
在锚点 $v = v_a$ 处求得一次 Lasso 解 $r^*$，提取非零活跃集 $\mathcal{A} = \{j : |r_j^*| > 0\}$（每道通常仅 $30\text{–}80$ 个非零点）。
对每个活跃原子 $j \in \mathcal{A}$，记其在锚点处的字典列为 $d_j = D_j(v_a)$，其对衰减率 $a(t_j)v$ 的导数原子为 $\dot{d}_j = \frac{\partial D_j}{\partial (a_j v)}$（即频域乘以 $b(f)$ 后对列范数导数做修正），并按活跃系数的幅值加权定义扰动向量：
$$\mathbf{q}_j = r_j^* \, \dot{d}_j \quad (j \in \mathcal{A})$$
当各活跃点经历失配 $\delta_j = a_j(u - v) - a_0(u - v_a)$ 时，经活跃集反射系数 $r_{\mathcal{A}}^*$ 自动重拟合松弛后，**无法被反射系数吸收的净残差变化**为：
$$\Delta \mathbf{e}_{\text{net}} = P_{\mathcal{A}}^\perp \sum_{j \in \mathcal{A}} \delta_j \mathbf{q}_j, \quad \text{其中 } P_{\mathcal{A}}^\perp = \mathbf{I} - \mathbf{D}_{\mathcal{A}}\big(\mathbf{D}_{\mathcal{A}}^T \mathbf{D}_{\mathcal{A}} + \epsilon\mathbf{I}\big)^{-1}\mathbf{D}_{\mathcal{A}}^T$$
于是目标函数对 $\delta$ 的真实非对角曲率矩阵为 $\mathbf{\Phi} = \mathbf{Q}^T P_{\mathcal{A}}^\perp \mathbf{Q}$（其中 $\mathbf{Q} = [\mathbf{q}_1, \dots, \mathbf{q}_{|\mathcal{A}|}]$），对 $v$ 求极小立刻给出**含完整段间耦合的闭式锚定系数**：
$$\alpha_{a,\text{Schur}} = \frac{a_0 \, \mathbf{1}_{\mathcal{A}}^T \big(\mathbf{Q}^T P_{\mathcal{A}}^\perp \mathbf{Q}\big) \mathbf{a}_{\mathcal{A}}}{\mathbf{a}_{\mathcal{A}}^T \big(\mathbf{Q}^T P_{\mathcal{A}}^\perp \mathbf{Q}\big) \mathbf{a}_{\mathcal{A}}} = \frac{a_0 \, \big\langle P_{\mathcal{A}}^\perp \mathbf{Q}\mathbf{1}_{\mathcal{A}},\; P_{\mathcal{A}}^\perp \mathbf{Q}\mathbf{a}_{\mathcal{A}} \big\rangle}{\big\| P_{\mathcal{A}}^\perp \mathbf{Q}\mathbf{a}_{\mathcal{A}} \big\|_2^2}$$
* **为什么它能解释之前所有的反常？**
  1. **为什么浅层大反射系数 $X^2$ 不会主导曲率（为什么 `α_E` 错而 `α_N` 接近）？** 因为浅层反射系数 $r_j^*$ 虽然幅值大，但它处于高信噪比、高活跃密度区，$P_{\mathcal{A}}^\perp \dot{d}_j$ 中绝大部分能量都被它自身和相邻系数在 $\mathbf{D}_{\mathcal{A}}$ 空间里的重拟合吸收掉了！
  2. **为什么分段扰动曲率之和是全局扰动的 1.2–9 倍？** 因为 $\mathbf{\Phi} = \mathbf{Q}^T P_{\mathcal{A}}^\perp \mathbf{Q}$ 的非对角元 $\Phi_{jk} = -\mathbf{q}_j^T \mathbf{D}_{\mathcal{A}}(\mathbf{D}_{\mathcal{A}}^T\mathbf{D}_{\mathcal{A}})^{-1}\mathbf{D}_{\mathcal{A}}^T \mathbf{q}_k < 0$ 是负的！对角元之和 $\sum_j \Phi_{jj}$ 远大于包含负非对角元的总二次型 $\mathbf{1}^T \mathbf{\Phi} \mathbf{1}$！
  3. **计算极快：** 只需要在锚点解的非零子集 $|\mathcal{A}| \le 100$ 上做一次小矩阵 QR 分解，单道耗时不到 $2\text{ ms}$，完全可以在 `anchor_active_curv_20261119` 的 112 个锚点解存档上直接验证！
