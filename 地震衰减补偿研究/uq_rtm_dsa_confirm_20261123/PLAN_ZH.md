# 预注册确认方案：倾角感知活排列镜面门控 UQ-RTM（`DSA-UQW1`）及浅层锚定实测修复验证（`uq_rtm_dsa_confirm_20261123`，开跑前提交）

## 0. 定位与防数据复用声明

本目录为 `uq_rtm_spg_20261123`（开发期 `G1–G4` 与 `G2` 事后机制核验）之后的**正式独立确认实验**。严格遵守：
- **零数据复用**：全部 4 套确认模型 `H1–H4`（`seed = 20261201–20261204`）与 1 套 Volve 式平盒实测复现模型 `H5_VolveAnc`（`seed = 20261205`）在构造形态、断层性质、低 Q 体分布、速度误差符号与随机种子上与 `G1–G4`、`D1–D4`、`C1–C8` 均完全独立、零重叠。
- **单一冻结候选方法**：本确认只冻结**单一候选算子 `DSA-UQW1`**（结构张量相干度加权倾角感知活排列镜面掩码＋倾角自洽双程 $t^*$，不除覆盖次数，不加易抖动的单频 Poynting 项），严禁在看分数后切换候选或调整任何参数/门槛。

## 1. 冻结算法定义（`DSA-UQW1`）

1. **第 1 遍 `R0` 与盲结构张量倾角估计（零额外波动方程求解）**：
   - 计算全炮未补偿叠加图像 $I_{R0}(\mathbf{x}) = \sum_f \operatorname{Re}\sum_s U_s^*(\mathbf{x}, f) V_s(\mathbf{x}, f)$，做常规 50 m 垂直平滑去除直流后求梯度 $(\partial_x I_{R0}, \partial_z I_{R0})$；
   - 高斯平滑（$\sigma = 4.0$ 格 $= 40\text{ m}$）构造二维结构张量 $J_{xx}, J_{zz}, J_{xz}$，计算原始倾角与相干度：
     $$\hat{\phi}(\mathbf{x}) = \operatorname{clip}\!\left(-\frac{1}{2}\operatorname{atan2}(2 J_{xz},\; J_{zz}-J_{xx}),\; -30^\circ,\; +30^\circ\right), \qquad c_{\text{st}}(\mathbf{x}) = \frac{\sqrt{(J_{zz}-J_{xx})^2 + 4 J_{xz}^2}}{\max(J_{zz}+J_{xx}, \epsilon)}$$
   - **相干度加权有效倾角** $\phi_{\text{eff}}(\mathbf{x}) = c_{\text{st}}(\mathbf{x}) \hat{\phi}(\mathbf{x})$（强反射层处 $c_{\text{st}}\approx 1$ 保持精确倾角，纯噪声/阴影区 $c_{\text{st}}\to 0$ 自动回退平层镜面，无自由调参）。
2. **倾角感知活排列镜面掩码 $w_{\text{DSA}}(\mathbf{x}, s)$ 与倾角自洽双程 $t^*$**：
   - 源入射角 $\theta_s = \operatorname{atan2}(x - x_s, \max(z - z_s, 10))$，反射出射角 $\theta_r = \operatorname{clip}(\theta_s + 2\phi_{\text{eff}}, -78^\circ, +78^\circ)$，地表镜面检波点 $x_{\text{mir}}(\phi_{\text{eff}}) = x + \max(z - z_s, 10)\tan\theta_r$；
   - 按各炮**真实活检波排列 $[x_{r,\min}(s), x_{r,\max}(s)]$** 计算入界距离 $d_{\text{in}} = \min(x_{\text{mir}} - x_{r,\min}(s), x_{r,\max}(s) - x_{\text{mir}})$，平滑余弦过渡带 `[-80, +120] m`，乘以近掠射抑制 `70°–80°`（不截断中大张角绕射/下穿射线），下限 `1e-4`；
   - 接收段 $t^*$ 亦按 $x_{\text{mir}}(\phi_{\text{eff}})$ 插值（与掩码几何自洽）。
3. **同域盲炮极性翻转估噪（`B = 整炮`）与期望算子 `UQW1`**：
   - 估噪 `iv` 在 $w_{\text{DSA}}$ 门控域内由炮翻转计算（天然兼容横向相关噪声）；
   - 贝叶斯参数与仓库完全一致：$\sigma = 0.447, \gamma = 1, s_{70} = \text{gain\_s}(70)$。
   - 同步计算严格公平对照：`RAW_{R0, IQ70, W1, UQW1}`、`SA_flat_{R0, IQ70, W1, UQW1}`（$\phi\equiv 0$ 消融）与 `DSA_{R0, IQ70, W1, UQW1}`。

## 2. 预注册确认模型矩阵（全部零反演犯罪，34 频点 `4:1:37 Hz`，5 档 $m\in\{0.5, 0.7, 1.0, 1.4, 2.0\}$）

### Part 1：4 套全新独立确认模型 `H1–H4`（`clean 15 dB` 与 `noisy 3 dB`）
1. **`H1`（向斜＋双偏心低 Q 透镜，`seed = 20261201`）**：向斜三界面＋双偏心低 Q 异常体（$Q=14$ 与 $Q=16$）。
2. **`H2`（逆断层＋大倾角反向陡倾层 `+22° / -16°`＋浅层低 Q 楔，`seed = 20261202`）**：比 `G2` 更陡（`22°`），严苛验证 `DSA-UQW1` 对陡倾构造与逆断层的无损增益及相对平层 `SA_flat` 的决定性修复。
3. **`H3`（Volve 式单边拖缆＋气烟囱＋12° 倾斜储层，`seed = 20261203`）**：非对称活拖缆 `[xs-860, xs+240]`（反向拖缆、~60% 排列外炮）叠加居中 $Q=15$ 气烟囱与 $12^\circ$ 倾斜目的层。
4. **`H4`（强横向相关噪声 `l30` $\sigma=60\text{ m}$＋`-3%` 反向速度误差＋横向 Q 梯度＋温和褶皱，`seed = 20261204`）**：与 `G4` 相反的速度偏低误差（`-3%`）叠加 60 m 横向相关噪声与横向渐变 $Q(x)\in[30, 110]$。

### Part 2：Volve 实测缺口修复验证模型 `H5_VolveAnc`（`seed = 20261205`）
模拟 Volve Stage 4 实测全链路：单边拖缆＋真实 Ricker 震源衰减正演，但偏移端**未知真实子波**而使用平盒 `W_mig = 1` 与 $(2\pi f)^2$ 成像条件，并在目的层抽 5 口虚拟井与**平盒带通真实垂向反射系数 `rw`** 做井震相关（严格复刻 `fld_eval.m`）。对比：
- `RAW_R0`（未补偿）
- `RAW_UQW1`（从地表 $z=0$ 补绝对 $t^*$，预期复现 Volve Stage 4 低于 `RAW_R0` 的缺口）
- **`DSA_UQW1_Anc`**（以浅层 `[190, 320] m` 窗口测得 `DSA_R0` 逐频振幅谱 $S_{\text{anc}}(f)$ 归一化去除 $(2\pi f)^2$ 与未知子波蓝倾，仅补偿相对 $\Delta t^* = \max(t^* - t^*_{\text{anc}}, 0)$）。

## 3. 预注册机械确认判据（`scripts/judge_dsa_confirm.py`，开跑前冻结）

主窗 `W_deep`（$x\in[180, 1420]\text{ m}, z\in[560, 950]\text{ m}$），参照为未门控全孔径弹性真值 `I_ref`：

- **`K1_DSA_vs_RAW`（全面超越现行 `RAW_UQW1`）**：在 `H1–H4` `noisy (3 dB)` `W_deep` 上，`DSA_UQW1 − RAW_UQW1` 满足：
  - 4 模型平均 $\Delta E_1 \ge +0.020$，且 **`4/4` 模型 $\ge +0.010$**；
  - 4 模型平均 $\Delta\text{SNR}_{\text{img}} \ge +1.00\text{ dB}$，且 **`4/4` 模型 $\ge 0.0\text{ dB}$**。
- **`K2_DSA_Fixes_SteepDip`（倾角感知在 22° 陡倾 `H2` 上严格胜过平层 `SA_flat` 且不伤 `RAW`）**：在 `H2` `noisy W_deep` 上，`DSA_UQW1 − SA_flat_UQW1` 满足 $\Delta E_1 \ge +0.008$、$\Delta E_{1,\text{HF}} \ge +0.010$、$\Delta\text{SNR}_{\text{img}} \ge +0.50\text{ dB}$；且 `DSA_UQW1 − RAW_UQW1` 在 `clean` 与 `noisy` 上均满足 $\Delta E_1 \ge 0.0$、$\Delta E_{1,\text{HF}} \ge -0.008$。
- **`K3_Fair_UQ_Robustness`（同门控公平贝叶斯稳健性）**：在 `H1–H4` `noisy W_deep` 上：
  - 过补偿 $m=0.5$ 上 `DSA_UQW1 − max(DSA_IQ70, DSA_W1)` 4 模型平均 $\ge +0.050$ 且 `4/4` 为正；
  - 相对全部同门控公平基线 `DSA_B+ = max(DSA_R0, DSA_IQ70, DSA_W1)`：`M2 (m=1.0) >= -0.025`（`4/4`），`M3 (m>1 worst) >= -0.035`（`4/4`）；
  - 最差情形 `min_m DSA_UQW1` 在 `4/4` 模型上同时严格高于 `min_m DSA_IQ70` 与 `min_m DSA_W1`。
- **`K4_Amp_NoHarm`（全模型振幅保真无损）**：在 `H1–H4` `noisy W_deep` 上，`DSA_UQW1` 的 `E_amp` 在 `4/4` 模型上 $\le \text{RAW\_UQW1} + 0.020$，且 4 模型均值 $\Delta E_{\text{amp}} \le 0.0$。
- **`K5_VolveAnc_WellTie`（Part 2 实测缺口修复检验，独立报告）**：在 `H5_VolveAnc`（`noisy`）的 5 井平均平盒井震相关上，`DSA_UQW1_Anc` 在全部 $m\in[0.5, 2.0]$ 均值与 $m=1.0$ 处均严格超过 `RAW_UQW1`（$\ge +0.05$）且不低于 `RAW_R0`（$\ge \text{RAW\_R0}$）。
