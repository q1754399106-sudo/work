# 方案：多炮匹配滤波贝叶斯 Q 补偿偏移（WS-UQ1 / WS-VA-UQ1）第 0 步核查与算子自检（T2.0，2026-11-23，开跑前提交）

## 0. 物理动机与问题定位

直接核查 `uq_rtm_20261117/scripts/migrate2.m:72` 与 `va_uq_field2_20261121/scripts/fe2_run.m:113` 发现：现行 `UQ-RTM` / `VA-UQ` 对每个成像点 $(z, x)$、频率 $f$ 上的 $n_s$ 炮采用**逐炮独立反演再等权累加**：
$$\hat{R}_{\text{old}} = \sum_{s=1}^{n_s} H_s^{\text{old}} C_s, \qquad H_s^{\text{old}} = \frac{A_{1,s} + s_{70} P_{1,s} C_{1,s} + P_{1,s}\,\text{iv}}{1 + s_{70} C_{1,s} + \text{iv}} \propto \exp(+\pi f t_s^*)$$
该架构在多炮异质衰减路径 $\{t_s^*\}$（如气烟囱绕射 vs 穿烟囱、Volve F-5 近道镜面 vs 孔径边缘大角度画弧）下违背多信道最大比合并（Maximal-Ratio Combining）原则：**对衰减最重、信噪比最低、画弧与速度失焦最严重的炮赋予了指数级最高权重**，并把叠后乐观信噪比 `iv` 套给单炮。

**新算子（`WS-UQ1` / `WS-VA-UQ1`）：**
令 $\bar{g}_s = \mathbb{E}_u[G(t_s^* u, f)]$（可由 `uq_tab` 的 $A_1/C_1$ 或 `va_tab` 直接查表得），炮归一化权重 $\sum_{s=1}^{n_s} a_s = 1$（主版 $a_s = 1/n_s$，对照源照明 $a_s = |U_s|^2/\sum|U_k|^2$），定义四个标量：
$$E_{\bar{g}} = \sum_s a_s |\bar{g}_s|^2, \quad L_{\bar{g}} = \sum_s a_s |\bar{g}_s|, \quad D_{\text{sig}} = \mathbb{E}_u\!\left|\sum_s a_s \bar{g}_s^* G_s(u)\right|^2, \quad D_{R0} = \mathbb{E}_u\!\left|\sum_s a_s G_s(u)\right|^2$$
则 $\gamma = 1$ 匹配滤波多炮贝叶斯算子为：
$$\hat{R}_{\text{WS-UQ1}}(z, x, f) = \frac{E_{\bar{g}} + (\text{iv}\cdot D_{R0} + s_{70}) \dfrac{E_{\bar{g}}}{L_{\bar{g}}}}{D_{\text{sig}} + (\text{iv}\cdot D_{R0} + s_{70}) E_{\bar{g}}} \sum_{s=1}^{n_s} \bar{g}_s^* C_s(z, x, f)$$
其满足三条严格边界恒等式：① $t_s^* \equiv 0 \implies \hat{R}_{\text{WS-UQ1}} \equiv R_0$；② 所有炮 $t_s^* \equiv t^*$ 相等时，逐位精确退化为 `UQW1`；③ 已知 $Q$ 极限下输出噪声方差相对逐炮反演严格缩减 $R_{\text{CS}} = \mathbb{E}_s[e^{-2\pi f t_s^*}]\mathbb{E}_s[e^{+2\pi f t_s^*}] \ge 1$（跳出单道定理 1 零和壁垒）。

## 1. 第 0 步两项零 LU 轻量核查（总内存 $< 500\text{ MB}$，不与同机 MATLAB LU 冲突）

1. **T2.0a（`scripts/step0_tstar_audit.m`）：真实程函炮间 $t_s^*$ 离散与定理 4 方差缩减比审计**
   - 只读读取已有缓存：
     - Volve 实测：`uq_rtm_20261117/cache/fld/ts.mat`（绝对 `TS`）与 `va_uq_field2_20261121/cache/anc/ts.mat`（2.27 km 锚定 `TS`）；
     - 合成开发集：`va_uq_step1_20261121/cache/fields_D2_vb.mat`（D2 气烟囱＋背斜）与 `fields_D4_vb.mat`（D4 浅层强衰减薄层）。
   - 在全部有效炮及物理孔径炮（$|x_s - x| \le 1.5 z$）两种口径下，统计：
     - 炮间标准差 $\operatorname{Std}_s(t_s^*)$、极差 $\max_s - \min_s$、分位数（P10/P50/P90）；
     - 定理 4 理论方差缩减比 $R_{\text{CS}}(f) = \mathbb{E}_s[e^{-2\pi f t_s^*}]\mathbb{E}_s[e^{+2\pi f t_s^*}]$（dB）在 $f = 15, 22, 35, 45\text{ Hz}$ 及 $m = 1, 0.5$ 下的中位数与 P90；
     - `UQW1_old` 与 `WS-UQ1` 对高 $t_s^*$ 炮（P75–P100）相对低 $t_s^*$ 炮（P0–P25）的相对权重比。
2. **T2.0b（`scripts/test_wsuq_identities.m`）：`ws_uq_apply.m` 算子三条边界恒等式与数值自检**
   - 验证：① $t_s^* = 0 \implies \max|\text{WS\_UQ1} - R_0| \le 10^{-6}$；② 均匀 $t_s^* = t^* \implies \max|\text{WS\_UQ1} - \text{UQW1}| \le 10^{-5}$；③ 异质 $t_s^*$ 下单频与宽带输出噪声方差、相关系数与独立 Python 校验一致。

## 2. 预登记门控（止损规则）

- **G0-1（算子数学自检）：** `test_wsuq_identities.m` 全部 3 条恒等式通过（误差 $\le 10^{-5}$）。
- **G0-2（真实炮间 $t_s^*$ 异质性门槛）：** 在活动物理孔径（$|x_s - x| \le 1.5 z$）内，以下两项至少满足一项才允许进入 T2.1 偏移运行：
  - (a) D2 气烟囱正下方窗（$x\in[950, 1450]\text{ m}, z\in[820, 1550]\text{ m}$）在 $f=35\text{ Hz}, m=1$ 下中位 $R_{\text{CS}} \ge 1.30$（$+1.1\text{ dB}$）；
  - (b) Volve F-5 储层窗（$z\in[2600, 3200]\text{ m}, x\in[-1.3, -1.1]\text{ km}$）在 $f=22\text{ Hz}, m=0.5$ 或 $f=35\text{ Hz}, m=1$ 下中位 $R_{\text{CS}} \ge 1.30$（$+1.1\text{ dB}$）。
- 若两项均 $< 1.30$，说明真实程函 $t_s^*$ 的炮间异质性太弱，立即停止 T2.1，不消耗偏移算力。
