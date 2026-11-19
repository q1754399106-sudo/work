# 叠前 Q 补偿的 AVO 误差预算与未动校道集炮检距谱比 Q 自校准

（《石油地球物理勘探》投稿初稿，2026-11-19。作者、单位、基金项目待填。图表见 `figures/` 与 `TABLES_ZH.md`；每个数字的出处见文末“数字出处”。）

---

**摘要：**叠前 Q 补偿常被用来提高分辨率并保持 AVO 特征，但 Q 模型本身的误差对 AVO 梯度的影响一直缺少定量分析。本文从局部平稳近似出发，推导了补偿后 AVO 振幅增益 Γ(x) 的无拟合参数闭式解，以及梯度截距比 G/A 对 Q 绝对水平的一阶灵敏度 ΔG/A ≈ −k(u − 1)，其中 u 为真实与模型 t* 之比，k 由通带有效频率、炮检距附加 t* 随 sin²θ 的增长率和反射层自身的 B/A 决定。由此得到 AVO 误差预算：补偿后逐事件 G/A 误差约为 k·σ_lnQ，补偿算子的设计消除不了这部分误差；σ_lnQ 超过约 0.2～0.45 时，最朴素的一步 Tikhonov 补偿已与最好的炮检距一致算子相当。对 Q 不确定度取期望的稳健补偿，因远道 t* 的绝对不确定度更大，会在 Q 正确时自行产生负的假 AVO 梯度（24～30 dB 时约 0.2～0.3）。为缩小 σ_lnQ，本文在未动校 CMP 道集上沿时距曲线取窗计算炮检距谱比，估计 Q 模型的全局校正系数 ĉ，并指出远道窗必须与近道窗等长：随时差压缩的远道窗会抹平谱衰减，产生随子波频带变窄而增大的正偏差（40 Hz 子波约 +0.08）。在 4 组从未使用过的模型上，修正后的估计器在 24 dB、合并 50 个 CMP 时 ln ĉ 的均方根误差为 0.091。端到端试验表明，用 ĉ 校准 Q 后，Q 偏离 2 倍以内时逐事件 G/A 误差由 0.4～1.8 降至 0.06～0.18；校准后各种沿炮检距路径的补偿算子在 AVO 上几乎没有差别，再次说明瓶颈在 Q 而不在算子。代价是模型 Q 本来正确时 G/A 误差增加约 0.06。全部结论来自合成数据。

**关键词：**叠前 Q 补偿；AVO；品质因子；炮检距谱比；误差预算；稳健补偿

**中图分类号：**P631　　**文献标识码：**A

---

## 0 引言

地层吸收使地震波高频成分随传播路径衰减。对叠前道集而言，同一反射点在不同炮检距上的射线路径长度不同，远道多经历的吸收会叠加在 AVO 响应上。只按垂向 t* 做时间域补偿会破坏 AVO，必须沿射线路径补偿，这一点已是共识[1-4]。围绕“沿路径补偿”，已有工作主要沿两个方向展开：一是改进补偿算子，如时间与炮检距两步反演[2]、射线路径 L1 正则化反演[3]、限频的数据驱动反 Q[5]、τ-p 域稀疏反演[6]和谱整形正则化[7]；二是改进叠前 Q 估计，如 QVO[8]、PSQI[9]和未动校 CMP 道集上的峰值频移法[10]。

两个方向通常各自评价：比较算子时假定 Q 已知，比较 Q 估计方法时不再追问估计误差传到 AVO 上有多大。已有研究对 Q 误差的处理停留在 ±15%～20% 的数值扰动对比[2,3]，没有给出 Q 误差与 AVO 梯度误差之间的定量关系。于是一个对实际处理很重要的问题没有答案：对 AVO 来说，瓶颈究竟是补偿算子，还是 Q 本身？

本文的工作包括四点：

1. 推导补偿后 AVO 振幅增益 Γ(x) 的闭式解及 G/A 对 Q 绝对水平的一阶灵敏度；
2. 由灵敏度给出 AVO 误差预算，并在同一框架下比较算子改进与 Q 精度的贡献；
3. 指出对 Q 不确定度取期望的稳健补偿会自行产生负的假 AVO 梯度；
4. 在未动校 CMP 道集上用炮检距谱比自校准 Q，分析远道窗被时差压缩造成的偏差并给出等长窗修正，再端到端检验校准对 AVO 的作用。

全部检验在合成数据上进行。所有判定规则在看结果之前写定，新模型只在一次检验中第一次使用；未通过的判定照实报告（表 3）。

## 1 方法原理

### 1.1 叠前道的路径吸收与补偿

层状介质中，零炮检距双程时间为 t 的反射，在炮检距 x 上经过的累积衰减为

$$a(t,x)=\int_{\text{ray}}\frac{\mathrm{d}s}{v\,Q}=a_0(t)+a_x(t,x) \tag{1}$$

式中 a₀ 为垂向双程 t*，a_x 为炮检距射线比垂向多出的部分。采用 Kolsky 型常 Q 频散模型[11]，路径传递函数为

$$G(f;a)=\exp\!\left[-\pi f a+\mathrm{i}\,2f\ln\!\left(f/f_{\mathrm r}\right)a\right] \tag{2}$$

参考频率 f_r 取 100 Hz。对每道（动校拉平后）记 G 为按 a(t, x) 逐时变化的衰减矩阵，广义 Tikhonov 补偿为

$$\left(\mathbf G^{\mathrm T}\mathbf G+\lambda\mathbf L^{\mathrm T}\mathbf L\right)\hat{\mathbf s}=\mathbf G^{\mathrm T}\mathbf d \tag{3}$$

L = I 时为一步 Tikhonov。L 取由炮检距附加吸收构造的零相位算子（符号为 exp[−π min(f, f_n) a_x]，f_n 由噪声放大上界确定）时为本文对比用的频率饱和阻尼，它属于广义 Tikhonov 的一种参数化，不作为本文的创新点。正则参数 λ 用残差原则[12]选取：在 log₁₀λ = 0, −0.5, …, −6 上由大到小搜索，取第一个使残差不超过噪声能量的 λ（噪声标准差已知）。

### 1.2 补偿后 AVO 增益的闭式解

设事件的真实路径吸收为 a_t，补偿使用模型路径 a_m。在局部平稳近似下，补偿后事件频谱乘以净传递

$$N(f,x)=\frac{G(f;a_{\mathrm t})\,G^{*}(f;a_{\mathrm m})}{|G(f;a_{\mathrm m})|^{2}+\lambda|L(f,x)|^{2}} \tag{4}$$

用子波振幅谱 S(f) 作匹配滤波读取振幅，事件在炮检距 x 上的振幅为 r(x)Γ(x)，r(x) 为反射系数，

$$\Gamma(x)=\frac{\int S^{2}(f)\,\mathrm{Re}\,N(f,x)\,\mathrm df}{\int S^{2}(f)\,\mathrm df} \tag{5}$$

Γ(x) 没有拟合参数。把 r(x)Γ(x) 按 Shuey 两项式[13]对 sin²θ 拟合得到 A′、B′，补偿引起的 G/A 偏差为 ΔG/A = B′/A′ − B/A。

### 1.3 G/A 对 Q 绝对水平的一阶灵敏度

记 u = a_t/a_m = Q_模型/Q_真（各层同比例），c = 1/u = Q_真/Q_模型。用炮检距一致的补偿算子时各道通带传递 ρ(f) = |G|²/(|G|² + λ|L|²) 相同，Q 有误时道 x 只多乘 exp[−πf(u−1)a_x(x)]，于是

$$\ln\Gamma(x)\approx \text{const}-\pi f_{\mathrm{eff}}(u-1)\,a_x(x),\qquad f_{\mathrm{eff}}=\frac{\int W\rho f\,\mathrm df}{\int W\rho\,\mathrm df} \tag{6}$$

W 为子波谱与带通的平方积。设 a_x ≈ η sin²θ，并以 sin⁴θ ≈ α + β sin²θ 在拟合的角度范围内线性化，代入 (A + B sin²θ)·exp(−πf_eff(u−1)η sin²θ) 并重新拟合，得

$$\Delta G/A\approx -k\,(u-1),\qquad k=\pi f_{\mathrm{eff}}\,\eta\left[1+\frac{B}{A}\beta-\left(\frac{B}{A}\right)^{2}\alpha\right] \tag{7}$$

方括号中的项来自拟合 (A + B sin²θ)Γ 时多出的 B sin⁴θ，**依赖反射层自身的 B/A**；略去它时，预测与数值的 R² 只有 0.17。式（7）的物理含义是：Q 的绝对水平每错 10%，逐事件 G/A 偏移约 0.1k；模型 Q 偏高（u > 1）时远道补偿不足，梯度偏负。

**误差预算。**若 ln u 的不确定度为 σ_lnQ，则在线性区，|ΔG/A| 的典型值（68% 分位数）约为 k·σ_lnQ。在本文模型上 k 为 0.4～1.9。这部分误差来自远道多出的 t* 误差 (u − 1)a_x，它与 AVO 梯度在单道上无法区分，任何不知道 u 的炮检距校正都消不掉。

### 1.4 对 Q 不确定度取期望的稳健补偿

叠后 Q 补偿中，可对 Q 的不确定度取期望构造稳健算子（期望 Wiener 形式，源自随机传递函数下的最小二乘滤波[14]）：

$$H=\frac{\mathrm E_u[G(f;u\,a_{\mathrm m})]^{*}}{\mathrm E_u|G(f;u\,a_{\mathrm m})|^{2}+\lambda},\qquad \ln u\sim N(0,\sigma^{2}) \tag{8}$$

推广到叠前道时，锚点取在地表，事件路径 t* 的绝对不确定度 ∝ a(x)，随炮检距增大。远道上 E|G|² 相对 |E G|² 增大得更快，算子补得更少，于是 Γ(x) 随炮检距下降，**即使 Q 给对（u = 1）也会产生负的假梯度**。记其大小为

$$D=\Delta G/A\big|_{\sigma}-\Delta G/A\big|_{\sigma=0} \tag{9}$$

同一 λ、同一子波、同一事件。D 只来自 σ，与 Q 是否正确无关。

### 1.5 未动校道集上的炮检距谱比 Q 自校准

AVO 与频率无关，同一反射层在近、远道上的反射系数谱形相同，因此近远道谱比可以消去反射系数谱色，只留下炮检距附加衰减。这正是叠前道集中关于 Q 绝对水平的信息。

为避免动校拉伸改变谱形，在未动校 CMP 道集上沿双曲时距曲线取窗。对零炮检距时间 t₀ 处长 T 的窗 [t₀ − T/2, t₀ + T/2]，近道功率谱为 P₀(f)；道 x 上取以 t_x(t₀) 为中心的窗，功率谱为 P_x(f)。两者都减去由噪声方差估计的噪声功率，并对 N 个 CMP 平均。于是

$$y=-\frac{\mathrm d}{\mathrm df}\left[\tfrac12\ln\frac{P_x(f)}{P_0(f)}\right]\approx \pi\,\Delta a_{\mathrm t}(t_0,x),\qquad \Delta a_{\mathrm t}=a_x^{\mathrm m}/c \tag{10}$$

斜率在自适应频带内拟合（12～90 Hz 内近、远道功率都超过噪声 3 倍、自低端连续的频段）。由模型 Q 计算 x = π⟨a_x^m⟩（窗内 Hann² 加权平均），对所有窗和道按最小二乘求全局系数：

$$\hat c=\frac{\sum x^{2}}{\sum x\,y} \tag{11}$$

时窗取 t₀ = 0.35～1.20 s、窗长 0.25 s、步长 0.1 s；拉伸比 t_x/t₀ > 1.3 的道不用。

**远道窗必须与近道等长。**若远道窗取为时距曲线上 [t_x(t₀ − T/2), t_x(t₀ + T/2)]，由 dt_x/dt₀ = t₀/t_x，窗长约为 T·t₀/t_x，随炮检距变短。Hann 窗越短，频谱平滑越宽，远道谱的高频衰减被抹平，Δa 被低估，ĉ 偏大；子波频带越窄，抹平的相对影响越大（2.4 节）。本文的估计器把远道窗取为以 t_x(t₀) 为中心、与近道等长的窗，并且不对功率谱做平滑。

校准后的 Q 模型为 ĉ·Q_模型，再接任何沿路径的补偿算子（图 1）。

## 2 模型试算

### 2.1 模型与检验程序

试算模型为 4 层介质（表 1），炮检距 0～1500 m、间隔 150 m，反射界面 t₀ = 0.38、0.65、0.90 s，零相位 Ricker 子波主频 60 Hz 和 40 Hz，反射系数用精确 Zoeppritz 方程计算。AVO 道集有孤立事件和薄层（每个事件后 16 ms 加反极性薄层）两类，拉平道集为判定场景，带双曲时差并动校（拉伸比 1.3 切除）的道集为保护场景。噪声为白噪声，按每道均方根定信噪比。真实 Q 等于 c 乘模型 Q，各层同比例。谱比估计用的超道集为 30% 稀疏反射系数，A ~ N(0, 1)，B = −0.8A + 0.6ε。

AVO 指标为有噪逐事件 G/A 误差：补偿结果带通到 10 Hz～子波主频后逐样点拟合 A、B，在每个事件 ±3/f_p 窗内求 B/A 与真值之差，三个事件取绝对值平均。另报告波形相对误差。

检验分 6 次进行，每次在运行前写定判定规则；用过的模型此后只作开发，新的检验使用从未用过的模型（表 1）。两处门槛是看过开发值后定的，在表 3 中注明。

**表 1  模型与检验**（详表见 `TABLES_ZH.md` 表 1a～1c）

| 项目 | 设置 |
|---|---|
| 层参数 | v_P 2100/2250/2450/2800 m/s，v_S 550/700/850/1000 m/s，ρ 1.8/2.1/2.3/2.4 g/cm³ |
| 观测系统 | 单 CMP，11 道，0～1500 m；采样 2 ms，记录 1.2 s |
| Q 模型 | 每层一个 Q，共 20 组（M12～M31），真 Q = c × 模型 Q，c = 0.5～2.0 |
| 子波 / 信噪比 | Ricker 60、40 Hz；12、18、24、30 dB |
| 补偿算子 | 频率饱和阻尼、分步、一步 Tikhonov、限频反 Q、只按垂向 t* 的 Tikhonov、期望算子（σ = 0.447） |
| 稳健补偿代价检验 | M12～M15，u = 1，λ = 10⁻⁴、10⁻³、10⁻² |
| 谱比第一次 / 第二次 | M16～M19 / M20～M23，c = 0.7、1.0、1.4，24 dB，N = 50 |
| 端到端 | M24～M27，c = 0.5～2.0，24、30 dB，10 个种子 |
| 校准后算子比较 | M28～M31，同上，6 个种子 |

### 2.2 闭式增益与灵敏度的检验

图 2 比较了由式（5）预测的逐事件 G/A 偏差与按时变算子数值计算的结果。一步 Tikhonov（σ = 0，Q 正确，72 点）回归斜率 0.999、R² 1.000；期望算子（σ = 0.447，u 取 0.5～2，504 点）斜率 0.974、R² 0.999。闭式对两类算子都成立，后文的误差预算可以解析进行。

式（7）的灵敏度在 17 个 Q 模型、2 种子波、3 个 λ、3 个事件上与闭式 Γ 的数值斜率对照：中、深事件（0.65、0.90 s）预测对数值的斜率 0.995、R² 0.997；浅事件（0.38 s，sin²θ 到 0.78）斜率 1.76、R² 0.74，推测是大角度远道 a_x 大、算子的噪声饱和使远道未被完全校正（未单独验证）。在独立正演、精确 Zoeppritz、残差原则选 λ 的端到端数据上（2.5 节，不校准的频率饱和阻尼，192 点），实测对预测的斜率 1.23、R² 0.75，符号全部正确；c = 0.5 时线性公式低估约 30%。因此式（7）可作为误差预算公式使用，精度约 ±30%，浅层大角度和强衰减时低估。

### 2.3 AVO 误差预算：算子消不掉 Q 的误差

图 3 给出先验 ln u ~ N(0, σ²) 下 |ΔG/A| 的 68% 分位数随 σ 的变化（λ = 10⁻³·⁵，闭式）。炮检距一致的算子（相当于频率饱和阻尼）在 σ = 0 时几乎无偏，误差随 σ 线性增长：40 Hz 子波 σ = 0.1、0.2、0.3 时约为 0.09、0.18、0.27，与 k·σ 一致。一步 Tikhonov 在 σ = 0 时自身已有 0.11～0.47 的偏差，但随 σ 增长较慢。两条曲线在 σ ≈ 0.2（40 Hz）至 0.45（60 Hz）处相交：**超过这一不确定度，最朴素的一步 Tikhonov 已与最好的炮检距一致算子相当。**σ = 0.447 时各算子误差为 0.27～0.57。这里不用均方根，因为 G/A 是比值，u 极端时截距趋零、比值发散，均方根被尾部主导。

开发中还试过两种在 Q 先验下“保 AVO”的算子（期望响应之比；零炮检距用期望算子、炮检距增量用确定性校正），σ = 0.447 时 68% 分位数比一步 Tikhonov 最多好 0.02，40 Hz 时差 0.08～0.13（开发数据，未作冻结检验）。这说明对 AVO 而言，算子设计的回报上限很低，能降低误差的只有缩小 σ_lnQ。

### 2.4 稳健补偿的 AVO 代价

在 M12～M15 上（Q 正确，u = 1，72 个单元），期望算子相对同 λ 的 Tikhonov 产生的假梯度 D 的中位数为 −0.161，94.4% 的单元为负；按 λ 分组，λ = 10⁻⁴、10⁻³、10⁻² 时分别为 −0.299、−0.163、−0.052；40 Hz 子波（−0.204）大于 60 Hz（−0.088）。Q 给对时，期望算子在 Tikhonov 原有偏差（0.39）之上再加约 0.16，总偏差达 0.56。在 u = 0.7～1.4 的范围内，期望算子的 AVO 误差都大于 Tikhonov，只在 u = 0.5、2 两端较小。

D 强烈依赖 λ，所以要看残差原则在各信噪比下实际选出的 λ（图 4）。每增加 6 dB，选出的 λ 约下降半个数量级。拉平道集上，D 的中位数在 30 dB 时为 0.29，24 dB 时为 0.21，18 dB 时为 0.13，12 dB 时为 0.06；动校道集上的 λ 大半档，D 相应更小。12～18 dB 时一步 Tikhonov 自身的 |ΔG/A| 为 0.42～0.75，比 D 大 5～20 倍。因此**稳健补偿的 AVO 代价在高信噪比道集（约 24 dB 以上，如超道集或偏移后道集）中才显著，量级 0.2～0.3**；噪声较大时被 Tikhonov 自身偏差淹没。

### 2.5 炮检距谱比的精度与等长窗修正

**第一次检验未通过。**在 M16～M19 上（24 dB，有 AVO，N = 50，768 个估计），原估计器（远道窗随时差压缩、功率谱 9 点平滑）的 ln(ĉ/c) 均方根为 0.109，最差格平均偏差 0.161，均超过预先写定的 0.10。失败几乎全来自 40 Hz 子波约 +0.08 的系统偏差，它在 30 dB、18 dB、无 AVO、N = 200 时都几乎不变，说明是估计器本身的偏差；60 Hz 子波基本无偏。AVO 的影响很小（无 AVO 时 0.102，有 AVO 时 0.109）。

**偏差来源。**无噪声、无 AVO 时（图 5），原估计器的 ln(ĉ/c) 在 30、40、60 Hz 子波下分别为 +0.353、+0.121、+0.033；远道窗改为等长并去掉平滑后为 +0.052、+0.011、0.000。改变拟合频带或炮检距附加 t* 的取法都无效。逐窗实测与预测的 Δt* 之比随炮检距下降（40 Hz：300 m 1.02，750 m 0.92，1200 m 0.85），与远道窗随炮检距变短一致。

**修正后在新模型上检验。**只改估计器（远道等长窗、不平滑），门槛不变，在从未使用过的 M20～M23 和新种子上重新检验：ln(ĉ/c) 均方根 0.091（平均 +0.021，标准差 0.088），最差格平均偏差 0.051，无失败估计，三项判定全部通过（图 6a）。描述性结果：30 dB 时 0.076，18 dB 时 0.108；N = 10 时 0.20（不够），N = 200 时 0.051（图 6b）；速度误差 ±2% 时 0.096、0.095；剩余偏差约 +0.02，原因未查。

按式（7）与图 3，σ_lnQ ≈ 0.09 对应逐事件 G/A 误差约 0.05～0.13。

### 2.6 端到端：校准对 AVO 的作用

在从未使用过的 M24～M27 上，c 取 0.5、0.7、1.0、1.4、2.0。每个单元用一个 50 个 CMP 的超道集估出 ĉ，再对同一 Q 模型的 AVO 道集用 ĉ·Q_模型 做频率饱和阻尼补偿，与三种参照比较：真值 Q（上限）、模型 Q（不校准）和模型 Q 的一步 Tikhonov。

**表 2′  端到端结果（拉平道集，单元中位数，每单元 10 个种子）**

| c | 真值 Q | 不校准 | 不校准 Tikhonov | **校准后** | 波形：不校准 / 校准后 | ln(ĉ/c) 均方根 |
|---|---|---|---|---|---|---|
| 0.5 | 0.041 | 1.758 | 1.764 | **0.176** | 0.915 / 0.563 | 0.096 |
| 0.7 | 0.023 | 0.639 | 0.713 | **0.108** | 0.703 / 0.432 | 0.069 |
| 1.0 | 0.014 | 0.014 | 0.080 | **0.080** | 0.300 / 0.330 | 0.069 |
| 1.4 | 0.010 | 0.398 | 0.349 | **0.065** | 1.022 / 0.282 | 0.087 |
| 2.0 | 0.008 | 0.582 | 0.565 | **0.061** | 2.909 / 0.242 | 0.125 |

Q 有误时，不校准的 G/A 误差为 0.4～1.8，校准后为 0.06～0.18，相对不校准降低 83%～90%（32 组的比值中位数 0.155）；波形误差也大幅降低，c = 2 时由 2.91 降到 0.24（不校准会过补偿、放大噪声）。按 5 个 c 平均，校准后 G/A 误差约 0.10，不校准约 0.68，真值 Q 为 0.019。ĉ 没有一次失败，c = 0.5 和 2.0 也能估得住。图 8 给出一次实现的 AVO 曲线示例。

**代价。**模型 Q 本来正确（c = 1）时，校准使 G/A 误差由 0.014 增加到 0.080，波形误差由 0.300 增加到 0.330，分别超过预先写定的 0.05 和 0.02，这两项判定未通过（表 3）。代价来自 ĉ 的估计噪声（c = 1 时 ln 均方根 0.069），按灵敏度推算是估计精度决定的下限，恰好等于用对 Q 的一步 Tikhonov（0.080）。

**动校道集。**带时差并动校的道集上，真值 Q、不校准、不校准 Tikhonov、校准后的 G/A 误差分别为 0.503、0.551、0.600、0.506：校准把 Q 引起的部分修到与真值 Q 相同，但动校拉伸带来约 0.5 的误差远大于 Q 的影响。

### 2.7 校准后的补偿算子比较

在从未使用过的 M28～M31 上（160 个单元，每单元 6 个种子），所有算子都使用同一个 ĉ·Q_模型（表 2）。沿炮检距路径的三种算子 G/A 误差几乎相同：频率饱和阻尼 0.097、分步 0.100、限频反 Q 0.099，配对差的中位数为 +0.0003 和 −0.0020，四分位都在 ±0.008 以内。一步 Tikhonov 为 0.123，劣势集中在 c ≤ 0.7（0.506、0.228）。只按垂向 t* 补偿为 1.095，比沿路径补偿大一个数量级，这一结果与已有共识[2,3]一致，在此作为基准对照。波形上限频反 Q 最好（0.245），分步最差（0.420）。

**表 2  校准 Q 后五种补偿算子的有噪逐事件 G/A 误差（单元中位数）**

| 算子 | 全部 | c = 0.5 | 0.7 | 1.0 | 1.4 | 2.0 | 波形误差 |
|---|---|---|---|---|---|---|---|
| 频率饱和阻尼 | 0.097 | 0.155 | 0.110 | 0.089 | 0.089 | 0.071 | 0.333 |
| 分步 | 0.100 | 0.178 | 0.122 | 0.090 | 0.086 | 0.070 | 0.420 |
| 限频反 Q | 0.099 | 0.180 | 0.121 | 0.091 | 0.088 | 0.071 | 0.245 |
| 一步 Tikhonov | 0.123 | 0.506 | 0.228 | 0.138 | 0.092 | 0.077 | 0.305 |
| 只按垂向 t* | 1.095 | 1.588 | 1.314 | 1.029 | 0.796 | 0.602 | 1.288 |

结合 2.3 节：Q 有误时任何算子都受 k·σ_lnQ 限制；Q 校准后，沿路径的算子之间又几乎没有差别。两头都说明，**对 AVO 而言，决定误差的是 Q 的绝对水平是否准确以及补偿是否沿炮检距路径，而不是补偿算子的具体形式。**

**表 3  预先写定的判定汇总**（详表与门槛来历见 `TABLES_ZH.md` 表 3）

| 检验 | 判定 | 结果 | 门槛 | 结论 |
|---|---|---|---|---|
| 稳健补偿代价 | \|D\| 中位数；负号比例；闭式斜率 | 0.161；94.4%；0.974 | ≥ 0.092；≥ 80%；0.9～1.1 | 通过 |
| 谱比，第一次 | 均方根；最差格；失败率 | 0.109；0.161；0 | ≤ 0.10；≤ 0.10；≤ 2% | **未通过**（前两项） |
| 谱比，第二次（新模型） | 同上 | 0.091；0.051；0 | 同上 | 通过 |
| 端到端 | 降低一半以上的组比例；距真值；Q 正确时 G/A 代价；波形代价；动校保护 | 100%；+0.070；+0.062；+0.027；92.2% | ≥ 90%；≤ 0.10；≤ 0.05*；≤ 0.02*；≥ 80% | **第 3、4 项未通过** |
| 校准后算子比较 | 与分步、限频反 Q 打平；垂向/沿路径；一步/沿路径 | +0.0003、−0.0020；10.7；1.35 | \|·\| ≤ 0.02*；≥ 5；≥ 1.2 | 通过 |

\* 看过开发值后定的门槛。

## 3 讨论

**Q 已正确时的代价。**校准在 c = 1 时增加约 0.06 的 G/A 误差，这是 ĉ 的估计噪声所致。ĉ 的离散随参与平均的 CMP 数下降（N = 200 时均方根 0.051），超道集更大时代价应按比例减小，但本文没有在端到端试验中检验。实际应用时是否校准，取决于对模型 Q 的把握：若模型 Q 来自可靠的 VSP 且误差在 10% 以内，校准的收益有限。

**动校拉伸。**动校道集上拉伸造成的 G/A 误差约 0.5，远大于 Q 的影响。所以在常规动校道集上做 AVO，主要问题不是 Q。本文的谱比估计在未动校道集上进行，不受拉伸影响；但 AVO 分析本身最好在偏移后的角道集上进行，这超出了本文范围。

**与已有方法的关系。**QVO[8]在动校后逐道调整窗长并校正拉伸，PSQI[9]转到 τ-p 域避开动校，峰值频移法[10]在未动校道集上沿反射轨迹取窗但未规定远近道窗长；张瑾等[15]的一维时窗试验表明，时窗短于子波时对数谱比法和中心频移法估出的 Q 远大于真值，但未涉及时差压缩。本文的估计器本身与这些方法同源，贡献在于指出未动校道集上远道窗随时差压缩造成的偏差机理及等长窗修正。Lazaratos 等[4,16]用由速度和炮检距算出的近远道相对 Q 补偿滤波（加动校拉伸补偿）均衡近远道频谱以保 AVO，Wu 和 Shi[17]面向 Q 补偿叠前时间偏移，在 CMP 道集上比较经过同一反射点的不同炮检距道，用谱比法和最小二乘求各炮检距的等效 Q，二者与本文“先估 Q、再沿路径补偿”的流程部分相近；本文的不同之处是等长窗修正和端到端的误差量化。限频反 Q[5]和广义 Tikhonov 反演[18]均有先例。

**稳健补偿的适用性。**叠后处理中，对 Q 不确定度取期望可以避免过补偿、提高稳健性；但到了叠前，同一思路会在 AVO 上付出代价。如果后续要做 AVO，叠前稳健补偿应慎用，或改为先缩小 Q 的不确定度。

**局限。**全部结论来自合成数据：模型浅（≤ 1.2 s）、双曲时差、各层 Q 同比例偏差（未测分层独立的形状误差）、超道集内 Q 无横向变化、残差原则使用已知噪声方差、无多次波和各向异性、子波不随炮检距变化。这些因素在实际资料上都会增加误差，实际应用前需要用 VSP 或井资料检验 ĉ。

## 4 结论

1. 补偿后 AVO 振幅增益有无拟合参数的闭式解，数值验证斜率 0.97～1.00。由此推出 G/A 对 Q 绝对水平的一阶灵敏度 ΔG/A ≈ −k(u − 1)，k 含反射层自身的 B/A，在本文模型上为 0.4～1.9。
2. AVO 误差预算：补偿后逐事件 G/A 误差约为 k·σ_lnQ，算子设计消除不了；σ_lnQ 超过约 0.2～0.45 时，一步 Tikhonov 与最好的炮检距一致算子相当。
3. 对 Q 不确定度取期望的叠前稳健补偿，在 Q 正确时自行产生负的假 AVO 梯度，24～30 dB 时约 0.2～0.3。
4. 在未动校道集上沿时距曲线取窗做炮检距谱比时，远道窗必须与近道等长；否则窗被时差压缩、谱被抹平，产生随子波频带变窄而增大的正偏差。修正后在 24 dB、50 个 CMP 时 ln Q 的误差约 0.09。
5. 用炮检距谱比校准 Q 后，Q 偏离 2 倍以内时逐事件 G/A 误差降低 83%～90%；校准后沿路径的各种补偿算子在 AVO 上几乎没有差别。代价是模型 Q 本来正确时 G/A 误差增加约 0.06。

## 参考文献

（按 GB/T 7714 排。2026-11-19 用网络检索核对过一轮，2026-11-20 在知网、本刊官网和 SEG 文献平台（GeoScienceWorld）原文页面核实了 [1][3][15][17]，结果见文末“参考文献核对记录”；仍标“待核”的条目投稿前必须在知网、SEG 文献平台或专利库中一手核实。排版时 [3] 题名中的“L1”的 1 用下标。）

[1] Luh P C. Wavelet dispersion and bright spot detection[C]//SEG Technical Program Expanded Abstracts 1988. Society of Exploration Geophysicists, 1988: 1217-1220. DOI: 10.1190/1.1892493.

[2] Li G, Liu Y, Zheng H, et al. Absorption decomposition and compensation via a two-step scheme[J]. Geophysics, 2015, 80(6): V145-V155.

[3] 程万里, 王守东, 孟巾钰, 王梓旭, 张俊杰. 基于L1范数正则化约束的叠前数据衰减补偿方法[J]. 石油地球物理勘探, 2023, 58(3): 567-579. DOI: 10.13810/j.cnki.issn.1000-7210.2023.03.008.

[4] Lazaratos S, Finn C. Deterministic spectral balancing for high-fidelity AVO[C]. SEG Technical Program Expanded Abstracts, 2004: 219-223.（题录取自他人文献的参考文献列表，原文待核）

[5] Ratnett N, Baldock S, Masoomzadeh H. A practical data-driven method for applying amplitude inverse Q that preserves AVO[C]. SEG Technical Program Expanded Abstracts, 2019: 620-624.

[6] Aharchaou M, Neumann E. Prestack Q compensation with sparse τ-p operators[J]. Geophysics, 2019, 84(5): V295-V305.

[7] Du Q, Wang W, Sun W, et al. Seismic attenuation compensation with spectral-shaping regularization[J]. Earth and Planetary Physics, 2022, 6(3): 259-274.

[8] Dasgupta R, Clark R A. Estimation of Q from surface seismic reflection data[J]. Geophysics, 1998, 63(6): 2120-2128.

[9] Reine C, Clark R, van der Baan M. Robust prestack Q-determination using surface seismic data: Part 1—Method and synthetic examples[J]. Geophysics, 2012, 77(1): R45-R56. DOI: 10.1190/GEO2011-0073.1.（副标题待核）

[10] Zhang C, Ulrych T J. Estimation of quality factors from CMP records[J]. Geophysics, 2002, 67(5): 1542-1547.

[11] Kolsky H. The propagation of stress pulses in viscoelastic solids[J]. Philosophical Magazine, 1956, 1(8): 693-710.

[12] Morozov V A. On the solution of functional equations by the method of regularization[J]. Doklady Akademii Nauk SSSR, 1966, 167(3): 510-512.（俄文原刊；英译本 Soviet Mathematics Doklady 的卷页未核实，故引原刊）

[13] Shuey R T. A simplification of the Zoeppritz equations[J]. Geophysics, 1985, 50(4): 609-614.

[14] Slepian D. Linear least-squares filtering of distorted images[J]. Journal of the Optical Society of America, 1967, 57(7): 918-922.

[15] 张瑾, 张国书, 王彦国, 李红星. 利用泰勒级数展开的振幅谱积分差值的Q值估计方法[J]. 石油地球物理勘探, 2022, 57(2): 320-330. DOI: 10.13810/j.cnki.issn.1000-7210.2022.02.009.

[16] ExxonMobil Upstream Research Company（发明人 Lazaratos S K）. Method for spectral balancing offset seismic data: WO2002065372A1[P]. 2002.（同族 US6516275B2、EP1360635A1；公开日与说明书内容待核）

[17] Wu J Z, Shi Y. Obtaining effective Q factors at different offsets from CMP records[J]. Geophysics, 2026, 91(4): S93-S100. DOI: 10.1190/GEO-2025-0353.

[18] 王本锋, 陈小宏, 李景叶, 陈增保. 基于反演的稳定高效衰减补偿方法[J]. 地球物理学报, 2014, 57(4): 1265-1274.（作者据英文版题录：Wang B F, Chen X H, Li J Y, Chen Z B. A stable and efficient attenuation compensation method based on inversion. Chinese Journal of Geophysics, 2015, DOI: 10.1002/cjg2.20181；中文版卷期页码待核）

---

## English abstract（初稿）

**AVO error budget of prestack Q compensation and Q self-calibration by offset spectral ratios on gathers without NMO correction**

**Abstract:** Prestack Q compensation is often applied to improve resolution while preserving AVO, but how errors in the Q model propagate into the AVO gradient has not been quantified. Under a local stationarity approximation we derive a closed-form, parameter-free AVO gain Γ(x) after compensation and the first-order sensitivity of the gradient-to-intercept ratio to the absolute Q level, ΔG/A ≈ −k(u − 1), where u is the ratio of true to model t* and k depends on the effective passband frequency, the growth of the offset t* with sin²θ, and the reflector's own B/A. This gives an AVO error budget: the per-event G/A error after compensation is about k·σ_lnQ and cannot be removed by operator design; for σ_lnQ above about 0.2–0.45 a plain one-step Tikhonov compensation performs as well as the best offset-consistent operator. A robust operator that takes the expectation over Q uncertainty creates a negative spurious AVO gradient even when Q is correct (about 0.2–0.3 at 24–30 dB), because the absolute t* uncertainty grows with offset. To reduce σ_lnQ we estimate a global Q correction from offset spectral ratios measured along the moveout curve on CMP gathers without NMO correction, and show that the far-offset window must have the same length as the near-offset window: a moveout-compressed far window smooths the spectral decay and causes a positive bias that grows as the wavelet band narrows (about +0.08 for a 40 Hz wavelet). On four unused models the corrected estimator gives an RMS error of 0.091 in ln ĉ at 24 dB with 50 pooled CMPs. End to end, calibration reduces the per-event G/A error from 0.4–1.8 to 0.06–0.18 when Q is wrong by up to a factor of two, and after calibration all path-consistent compensation operators give nearly the same AVO, confirming that Q, not the operator, is the bottleneck. The cost is an increase of about 0.06 when the model Q is already correct. All results are synthetic.

**Keywords:** prestack Q compensation; AVO; quality factor; offset spectral ratio; error budget; robust compensation

---

## 投稿前待办（不属于正文）

1. **图件：**~~改中文标签~~ **已做**（2026-11-19，`figures/`，英文版移到 `figures/en/`）。终稿字体（宋体 + Times New Roman）需在 Windows 上运行 `python scripts/make_figs.py zh` 重出一次；当前文件用的是度量兼容的 Liberation Serif + Noto Serif CJK SC。图题按本稿：图 1 流程，图 2 闭式 Γ 检验，图 3 误差随 σ，图 4 D 随 λ，图 5 远道窗压缩，图 6 ĉ 精度，图 7 端到端，图 8 示例道集。正文中图 7 的数字以表 2′ 给出，可二选一以控制篇幅。
2. **参考文献：**2026-11-19 已检索核对一轮（见“参考文献核对记录”）。[4][9][12][16][18] 已补全或改引；2026-11-20 已在原文页面核实 [1][3][15][17] 并补全题录（[1] 改引 Luh 1988 SEG 扩展摘要）。仍待办：[1] 正文只有 PDF（需订阅），引言“共识[1-4]”中 [1] 的具体内容未读到；[17] 只读到摘要，是否实际用分炮检距 Q 做了偏移补偿需读全文。
3. **篇幅：**本刊常见篇幅 8～12 页。若超长，可把 2.4 节的 λ 实测并入一句话，把表 3 移入附录。
4. **补充材料：**代码与预先写定的方案（各 PLAN 文件及提交号）可作为补充材料或数据可用性声明。
5. **作者、单位、基金、通讯作者信息**待填。

## 参考文献核对记录（2026-11-19，网络检索）

本会话只能用检索结果摘要核对，不能打开原文页面；以下“已核”指题录信息在至少一个独立来源中出现。

| 条目 | 结果 | 依据 |
|---|---|---|
| [1] Luh 1993 | **未找到** | 只找到 Luh 1989 AAPG 摘要 “Wavelet dispersion and bright-spot detection”；1993 年文献在他人综述中被引用，但题名、出处未见 |
| [3] 程万里等 2023 | **未找到** | 题名与查新报告给出的 DOI 均检索不到 |
| [4] Lazaratos & Finn 2004 | 已核（二次来源） | 2007 年 CSEG 摘要的参考文献列表：74th SEG Annual Meeting, Expanded Abstracts, 219–223 |
| [9] Reine 等 2012 | 已核 | 卷期页码与 DOI 一致；第 2 部分副标题 “3D case study” 已见，第 1 部分副标题未见原文 |
| [12] Morozov 1966 | 已核（改引原刊） | Math-Net.Ru 等：Doklady Akademii Nauk SSSR, 167(3): 510–512 |
| [15] 张瑾等 2022 | **未找到** | 只找到同一作者 2023、2024 年的相关论文 |
| [16] WO2002065372A1 | 已核（部分） | Google Patents 记录存在；ExxonMobil，发明人 Lazaratos；同族 US6516275B2、EP1360635A1；方法为近远道“相对 Q 补偿滤波”加拉伸补偿，是否由谱比估计 Q 未见 |
| [17] Geophysics 91(4):S93 | **未找到** | 无法确认该文存在；查不到就删除 |
| [18] 王本锋等 2014 | 已核（英文版） | Chinese Journal of Geophysics 英文版 2015，DOI 10.1002/cjg2.20181，作者 Wang, Chen, Li, Chen |
| 其余（[2][5]–[8][10][11][13][14]） | 未在本轮复核 | [2][5][7][8][10][14] 的卷期页码与查新报告一致；[6][11][13] 为常用文献，仍建议投稿前统一用 Crossref 核对 |

**2026-11-20 原文页面核实（取代上表中 [1][3][15][17] 的“未找到”）：**

| 条目 | 结果 | 原文页面与说明 |
|---|---|---|
| [1] Luh | **改引 1988** | 没有 Luh 单独署名的 1993 年文献。1993 年 Castagna & Backus 编《Offset-Dependent Reflectivity》（SEG）第 3 章 “Processing” 由 Resnick、Luh、Shang、McDonald、Gardner 合著（DOI 10.1190/1.9781560802624.ch3），官方页面不列章内分节，Luh 部分的内容无法确认。改引与 AAPG 1989 摘要同题的 SEG 1988 扩展摘要：GeoScienceWorld 页面 pubs.geoscienceworld.org/segeab/proceedings-abstract/SEGEAB.7/1/1217/672070，pp. 1217–1220，属 AVO 理论与分析专场；正文仅 PDF，需订阅，未读 |
| [3] 程万里等 2023 | 已核（知网） | link.cnki.net/doi/10.13810/j.cnki.issn.1000-7210.2023.03.008；作者 5 人，题录与原稿一致；摘要含 AVA 特征恢复、与叠后补偿和常规叠前反 Q 对比、Q 值敏感度分析 |
| [15] 张瑾等 2022 | 已核（本刊官网、知网、全文） | ogp-cn.com.cn 官网页 320–330（知网记为 320-330+243-244，后者应为彩页）；全文 html.rhhz.net/SYDQWLKT/2022-2-320.htm。2.1.2 节“时窗宽度试验”：一维合成记录、Q=50、10–100 Hz，时窗 10/20 ms（非完整子波）时 LSR、CFS 估出的 Q 远大于真值，LSAD、ASID 接近真值；正偏差限于 LSR、CFS，原文用“非完整子波”而非“截断”，正文表述已据此改写。2023, 58(6): 1423-1435 是同一作者的另一篇（DOI …2023.06.014），与本条无关 |
| [17] Wu & Shi 2026 | 已核（GeoScienceWorld 摘要页） | pubs.geoscienceworld.org/seg/geophysics/article-abstract/91/4/S93/730993；2026-06-24 在线，91(4): S93–S100。摘要：为 QPSTM 定义非零入射角的等效 Q（弯曲射线走时），比较 CMP 道集中经同一反射点的不同炮检距道，谱比法建方程、最小二乘求各炮检距等效 Q，合成与实际资料验证。摘要未明说用该 Q 做了偏移补偿，正文表述已改为“面向 QPSTM”；全文需订阅，未读 |

**风险提示：**[3][15][17] 已确认存在，原“查不到就删除”的风险解除。[1] 的 1988 扩展摘要只核实了题录，若审稿人质疑其内容与“沿路径补偿”共识的对应关系，可从 [1-4] 中去掉 [1]，由 [2]–[4] 支撑。

## 数字出处

| 位置 | 文件 |
|---|---|
| 2.2 闭式检验 | `uq_avo_20261118/RESULTS_GATE1_ZH.md`（P1）；灵敏度 `uq_avo_20261118/THEORY_K_ZH.md` |
| 2.3 误差预算 | `uq_avo_20261118/DEV_GATE2_ZH.md`（末节 68% 分位数更正） |
| 2.4 稳健补偿代价 | `uq_avo_20261118/RESULTS_GATE1_ZH.md`、`RESULTS_LAMHI_ZH.md` |
| 2.5 谱比 | `RESULTS_GATE_A1_ZH.md`、`DIAG_A1_ZH.md`、`RESULTS_GATE_A1B_ZH.md` |
| 2.6 端到端 | `RESULTS_GATE_A2_ZH.md`、`results/a2_verdict.json` |
| 2.7 算子比较 | `RESULTS_GATE_A3_ZH.md`、`results/a3_verdict.json`、`TABLES_ZH.md` 表 2 |
| 表 3 | `TABLES_ZH.md` 表 3 |
