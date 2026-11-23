# 预注册实验与工程合流方案：UQ-RTM 最终推荐统一流程 `DSA-UQW1b-Anc`（`uq_rtm_unified_20261124`，开跑前冻结）

- **日期**：2026-11-24
- **工作树与分支**：`D:\ChatGPT\研究\.claude\worktrees\uq-rtm-unified-20261124`（分支 `agy/uq-rtm-unified-20261124`，由 `origin/main` 独立切出，绝不污染任何已有工作树）
- **核心目标**：将两支已分别取得突破的技术主线——**`UQW1b`（$\gamma=1$ 期望成像条件 + $B=300$ 整炮极性翻转估噪）**与 **`DSA-UQW1-Anc`（盲结构张量倾角感知活排列镜面门控 + 参数化浅层锚定谱均衡 `PSAE` + 频带内速度一致参考频率 $f_{\text{ref}}=25\text{ Hz}$）**——合流为单一、自洽、可复用的最终推荐流程 **`DSA-UQW1b-Anc`**，并在**合成多物理联合盲测（`U1–U5`）**、**Volve 北海实测 4–25 Hz 全因子消融（`fld`）**与**Volve 北海实测 4–45 Hz 宽频 Stage 4H 合流（`fld45`）**三层体系上完成端到端预注册闭环验证。
- **只读与盲测铁律**：
  1. 对上游实测缓存 `uq_rtm_20261117/cache/fld`、`va_uq_field2_20261121/cache/inv300.mat` 与 `study-research-docs-plan-c4d4c6/.../uq_rtm_20261117/cache/fld45` 严格**只读**，绝不修改任何源文件；
  2. 算法全部参数（盲倾角 $\phi_{\text{eff}}$、镜面掩码 $w_{\text{DSA}}$、$B=300$ 估噪 $\text{iv}_{300}$、浅层锚定窗 `[1000, 1600] m` 斜率 $\beta_{\text{anc}}$ 与相对衰减 $\Delta t^*_{\text{rel}}$、$f_{\text{ref}}=25\text{ Hz}$）**100% 仅由地震数据、采集几何与偏移速度/$Q$ 模型决定，算法全程严禁触碰测井曲线 `F-5` / `F-14`**；
  3. 重型中间频率切片/波场缓存统一写入工作区外 `D:\tmp_uqrtm_unified_20261124\`，严禁提交至 Git 仓库。

---

## 1. 统一推荐流程 `DSA-UQW1b-Anc` 完整数学定义（`scripts/uq_rtm_unified_core.m`）

全流程零额外波动方程 LU 分解，由四大物理闭环模块串联组成：

### 模块 A：盲结构张量倾角提取与活排列防假频镜面门控（`DSA`）
1. **盲倾角估计**：对第 1 遍未补偿全孔径像 $I_{\text{RAW}, R_0}(z, x)$ 做 40–50 m 垂直高斯高通去直流后，计算梯度 $(\partial_x I, \partial_z I)$ 与高斯平滑二维结构张量 $\mathbf{J} = G_{\sigma_s} \ast (\nabla I \nabla I^T)$（实测 $\sigma_s = 60\text{ m}$，合成 $\sigma_s = 40\text{ m}$），闭式提取相干度加权有效倾角：
   $$\hat{\phi}(\mathbf{x}) = \operatorname{clip}\!\left(-\frac{1}{2}\operatorname{atan2}(2J_{xz},\; J_{zz}-J_{xx}),\; -30^\circ,\; +30^\circ\right), \quad c_{\text{st}}(\mathbf{x}) = \frac{\sqrt{(J_{zz}-J_{xx})^2 + 4J_{xz}^2}}{\max(J_{zz}+J_{xx},\, \epsilon)}, \quad \phi_{\text{eff}}(\mathbf{x}) = c_{\text{st}}(\mathbf{x})\hat{\phi}(\mathbf{x})$$
2. **倾角感知活排列镜面掩码 $w_{\text{DSA}}(\mathbf{x}, s)$**：
   对成像点 $\mathbf{x}=(z,x)$ 与第 $s$ 炮（炮点 $(z_s, x_s)$，该炮实际活检波器区间 $[x_{r,\min}(s), x_{r,\max}(s)]$）：
   $$\theta_s = \operatorname{atan2}(x - x_s, \max(z-z_s, 10)), \quad \theta_r = \operatorname{clip}(\theta_s + 2\phi_{\text{eff}}, -78^\circ, +78^\circ), \quad x_{\text{mir}} = x + \max(z-z_r, 10)\tan\theta_r$$
   计算入界距离 $d_{\text{in}} = \min(x_{\text{mir}} - x_{r,\min}(s), x_{r,\max}(s) - x_{\text{mir}})$，施加活排列余弦过渡与防炮采样假频角限（实测 OBC 50 m 炮距：$d_{\text{in}}\in[-150, 250]\text{ m}, \theta_{\max}=50^\circ, \Delta\theta=10^\circ$；合成小尺度模型：$d_{\text{in}}\in[-80, 120]\text{ m}, \theta_{\max}=75^\circ, \Delta\theta=10^\circ$，下限 $10^{-4}$，不做粗糙逐点覆盖除法）：
   $$w_{\text{DSA}}(\mathbf{x}, s) = \max\!\Big(T_{\cos}(d_{\text{in}})\, T_{\cos}(\theta_{\max}-|\theta_s|)\, T_{\cos}(\theta_{\max}-|\theta_r|),\; 10^{-4}\Big)$$
3. **倾角自洽双程衰减时间 $t^*_{\text{DSA}}(\mathbf{x}, s)$**：直接复用单程程函 $t^*$ 场 $F(\mathbf{x}, k)$，按倾角投影镜面点 $x_{\text{mir}}(\phi_{\text{eff}})$ 线性插值接收支 $t^*$ 与炮支相加（零额外程函求解）。

### 模块 B：整炮分块随机极性翻转盲估噪（`B = 300`，定理 3 推演至 DSA 门控域）
为每炮分配一个独立 Rademacher 随机符号 $\epsilon_s \in \{\pm 1\}$（同炮全部检波道共享 $\epsilon_s$，即块长 $B=300$）。由定理 3，整炮翻转使块内横向相关噪声保留率 $r_{300} = t$（彻底消除逐道 $B=1$ 在横向相关噪下的 13–28 倍低估），同时多炮相干反射信号在炮间随机翻转下相互抵消。在 $\pm 5\text{ Hz}$ 频带内求和并做空间方框平滑（实测 400 m，合成 160 m）：
$$p(\mathbf{x}, f_0) = \operatorname{box}_L\!\left(\Big[\sum_{|f-f_0|\le 5\text{ Hz}} I(\mathbf{x}, f)\Big]^2\right), \quad p_f(\mathbf{x}, f_0) = \operatorname{box}_L\!\left(\Big[\sum_{|f-f_0|\le 5\text{ Hz}} I_{\text{flip}}(\mathbf{x}, f)\Big]^2\right), \quad \text{iv}_{300} = \min\!\left(\frac{p_f}{\max(p-p_f, \epsilon)},\, 10^3\right)$$
本实验同时评估 **DSA 域内 `B=300`（`iv3_dsa`，在 $w_{\text{DSA}}$ 门控后执行带内求和翻转）** 与 **全孔径 `B=300`（`iv3_raw`，与 Stage 4H 严格同口径）**。

### 模块 C：参数化平滑浅层锚定谱均衡（`PSAE`）与相对衰减 $\Delta t^*_{\text{rel}}$
针对未知震源子波并采用平盒源 $W_{\text{mig}}=1$ 与 $(2\pi f)^2$ 成像条件的实测工作流（及合成实测孪生体 `U4`）：
1. 在浅层参考窗 $\Omega_{\text{anc}}$（Volve 实测 `[1000, 1600] m, x ∈ [-2500, 1500] m`；1 km 合成 `U4` `[190, 320] m`）提取锚定衰减基线与相对衰减：
   $$t^*_{\text{anc}} = \operatorname{median}_{\Omega_{\text{anc}},\, s}\big(t^*_{\text{DSA}}\big), \qquad \Delta t^*_{\text{rel}}(\mathbf{x}, s) = \max\!\big(t^*_{\text{DSA}}(\mathbf{x}, s) - t^*_{\text{anc}},\; 0\big)$$
2. 在 $\Omega_{\text{anc}}$ 内测量盲 `DSA_R0(f)` 的逐频 RMS 振幅谱 $S_{\text{anc}}(f)$，在有效信号频带 $[f_{\text{fit},1}, f_{\text{fit},2}]$（4–25 Hz 实测 `[6, 22] Hz`；4–45 Hz 实测 `[6, 38] Hz`；合成 `U4` `[6, 32] Hz`）做单参数对数线性拟合 $\ln S_{\text{anc}}(f) \approx a_0 + \beta_{\text{anc}} f$，构造无凹口平滑均衡权重：
   $$W_{\text{eq}}(f) = \exp\!\big(-\beta_{\text{anc}}(f - f_{\text{mid}})\big)$$
   *(对于震源子波完全已知且以无衰减弹性 RTM 像 $I_{\text{ref}}$ 为真值对照的常规合成测试 `U1–U3, U5`，启用 `anchor_mode='none'`，即 $W_{\text{eq}}\equiv 1$、直接补偿绝对 $t^*_{\text{DSA}}$；对于未知子波实测流与 `U4`，默认启用 `anchor_mode='psae'`。)*

### 模块 D：频带内速度一致（$f_{\text{ref}}=25\text{ Hz}$）$\gamma=1$ 贝叶斯期望成像条件（`UQW1b`）
在 $\hat{x} = \Delta t^*_{\text{rel}}(\mathbf{x}, s)/m$ 处，采用 41 点 Gauss–Hermite 对数正态求积（$\sigma = 0.447$，**$f_{\text{ref}} = 25\text{ Hz}$**，$s_{70} = \text{gain\_s}(70)$）计算 $A_1, C_1, P_1$，构造噪声主导区自动回退纯相位频散校正（$\lim_{\text{iv}\to\infty} H_1 = P_1$）的统一算子：
$$H_{\text{UQW1b}}(\mathbf{x}, s, f) = \frac{A_1 + s_{70}\, P_1\, C_1 + P_1\, \text{iv}_{300}(\mathbf{x}, f)}{1 + s_{70}\, C_1 + \text{iv}_{300}(\mathbf{x}, f)}$$
$$I_{\text{DSA-UQW1b-Anc}}(\mathbf{x}) = \sum_f W_{\text{eq}}(f)\, \operatorname{Re}\sum_{s=1}^{N_s} w_{\text{DSA}}(\mathbf{x}, s)\, H_{\text{UQW1b}}(\mathbf{x}, s, f)\, C_s(\mathbf{x}, f)$$

---

## 2. 三层全因子闭环验证设计（零数据复用、零反演犯罪）

### Tier 1：全新多物理联合合成盲测套件 `U1–U5`（`scripts/run_unified_synth.m`）
全 2D 频域黏声 Helmholtz 精确正演（`4:1:37 Hz`，34 频，100×160 网格，49 炮 × 160 检波器，$m \in \{0.5, 0.7, 1.0, 1.4, 2.0\}$），**全部采用实测频带内速度约定 $v_{\text{band}}(25\text{ Hz})$** 与**全新独立随机种子 `20261301–20261305`**，首次将四大难点同时耦合：
- **`U1`（气烟囱 + 倾角反射层 + 横向相关噪 `l30` $\sigma_r=50\text{ m}$ + $v_{\text{band}}(25\text{ Hz})$）**：同场检验 `B=1` vs `B=300`（RAW 与 DSA 域对比真值 `iv_true`）、`RAW` vs `DSA`、`f_ref=100` vs `25 Hz`；
- **`U2`（逆断层 + `+22°/-16°` 陡倾构造 + `l30` 相关噪 + $v_{\text{band}}(25\text{ Hz})$）**：检验陡倾下 `DSA` 相对 `SA_flat` 与 `RAW` 的无损增益及与 `B=300` 的协同；
- **`U3`（Volve 式单边拖缆 ~65% 排列外炮 + 气烟囱 + `l30` 相关噪 + $v_{\text{band}}(25\text{ Hz})$）**：检验高比例排列外画弧剔除与 `B=300` 的协同；
- **`U4`（Volve 实测全链路孪生体：未知 Ricker 子波 $\to$ 平盒 `W_mig=1` + $(2\pi f)^2$ + 单边拖缆 + `l30` 相关噪 + $v_{\text{band}}(25\text{ Hz})$ + 5 口精细薄层虚拟井）**：同时激活全部四大组件（`B=300` + `DSA` + `PSAE` + `f_ref=25 Hz`），按 `fld_eval` 同口径评价 5 井固定带 `[4–18 cyc/km]` 相关与储层谱斜率；
- **`U5`（高信噪比独立噪 + `-3%` 速度误差 + 横向 Q 梯度安全边界检验）**：验证高信噪比与速度失配下合流流程无害性。

### Tier 2：Volve 北海实测 OBC 4–25 Hz 全因子合流（`scripts/run_unified_volve25.m`）
严格只读 `cache/fld`（10 m 网格，85 频 `4.00:0.25:25.00 Hz`，124 炮）与 `va_uq_field2_20261121/cache/inv300.mat`。单遍求解 Helmholtz 并临时流式缓存 `Cs` 至 `D:\tmp_uqrtm_unified_20261124\cs_fld25`（运行结束即自动清理），完整计算 **DSA 域内 `B=300`（`iv3_dsa`）**、**全孔径 `B=300`（`iv3_raw`）**与 **`B=1`（`iv1_raw`）**，同场产出完整 2×2×2×3 因子对照（$m \in \{0.5, 0.7, 1.0, 1.4, 2.0\}$）：
- 基线与复现：`RAW_R0`（校验误差 $<10^{-4}$）、`DSA_R0`、`DSA_R0_Eq`、`RAW_UQW1_Unanc_25`、`RAW_UQW1b_Unanc_25`、`DSA_UQW1b_Unanc_25`、`RAW_UQW1b_PAnc_25`；
- 公平竞争者：`DSA_IQ70_PAnc_25`、`DSA_W1b_PAnc_25`（用 `iv3_dsa`）、`DSA_PH_PAnc_25`；
- 统一推荐主方法与相位/估噪消融：
  - **`DSA_UQW1b_PAnc_25`**（**主推**：`DSA` + `iv3_dsa` ($B=300$) + `PSAE` + $f_{\text{ref}}=25\text{ Hz}$）；
  - `DSA_UQW1bRaw_PAnc_25`（`iv3_raw` 对照）、`DSA_UQW1_PAnc_25`（`B=1` 对照）；
  - `DSA_UQW1b_PAnc_Amp`（纯振幅 $f_{\text{ref}}=f$ 伴随变体）、`DSA_UQW1b_PAnc_100`（遗留 $f_{\text{ref}}=100\text{ Hz}$ 伴随变体）。

### Tier 3：Volve 北海实测 OBC 4–45 Hz 宽频 Stage 4H 合流（`scripts/run_unified_volve45.m`）
严格只读 `study-research-docs-plan-c4d4c6/.../cache/fld45`（6.25 m 网格，165 频 `4.00:0.25:45.00 Hz`，124 炮，直接复用已缓存的 `inv.mat` [`iv` $B=1$ & `iv3` $B=300$]、`tfields.mat`、`ts.mat` 与 `images.mat: R0` 做盲倾角），单遍逐频可续跑，轻量化 2D 切片存 `D:\tmp_uqrtm_unified_20261124\slices_fld45`，在 Exact Stage 4H `m = [0.6, 0.85, 1.0, 1.2, 1.4]` 与固定评价带 `tie_fix`（8–40 Hz，`4.17–20.86 cyc/km`）、`tie_r0`、零时移相关 `L=0`、`F-14`、上覆画弧占比 `Smile_Frac`（`[1600, 2400] m`）及储层谱斜率 `Res_Slope` 下同场对比：
- 复现校验：`RAW_R0`（对标 Stage 4H `0.4740`）、`RAW_UQW1b_Unanc_100`（对标 Stage 4H `UQW1b = 0.5380`）；
- 消融与公平基线：`DSA_R0`、`DSA_R0_Eq`、`DSA_UQW1b_Unanc_25`、`DSA_UQW1b_Unanc_100`、`RAW_UQW1b_PAnc_25`、`DSA_UQW1_PAnc_25`（$B=1$）、`DSA_IQ70_PAnc_25`、`DSA_W1b_PAnc_25`；
- **主推合流与相位伴随**：**`DSA_UQW1b_PAnc_25`**（**主推**：`DSA` + $B=300$ + `PSAE` + $f_{\text{ref}}=25\text{ Hz}$）、`DSA_UQW1b_PAnc_Amp`（纯振幅）、`DSA_UQW1b_PAnc_100`（$f_{\text{ref}}=100\text{ Hz}$）及 `EnvAnc`（8 Hz 宽窗平滑浅层包络归一化）后处理对照；并同步累计 `I_dsa_f` 与 `I_dsa_flip300_f` 以定量审计 4–45 Hz 下 `DSA` 域内 vs `RAW` 域内 $B=300$ 估噪差异。

---

## 3. 预注册自动判决门槛（`scripts/judge_unified.py`，运行前冻结）

- **`G1_Theorem3_B300_In_DSA`**：在 `U1–U3` 横向相关噪 `l30` 下，`RAW B=1` 的 $\log_2(\text{iv}_{\text{est}}/\text{iv}_{\text{true}})$ 中位数 $\le -2.0$（严重低估 $\ge 4\times$），而 **`DSA B=300`** 在 `3/3` 模型上均满足 $|\log_2(\text{iv}_{\text{est}}/\text{iv}_{\text{true}})| \le 1.0$。
- **`G2_Synth_Confluence_Beats_Components`**：在 `U1–U3`（`noisy l30`, `v_band=25 Hz`, `W_deep`）上，统一 **`DSA_UQW1b_25`** 满足：
  1. 相对 `RAW_UQW1_25`（$B=1$, 无 DSA）：3 模均值 $\Delta E_1 \ge +0.025$ 且 `3/3 > 0`；
  2. 相对 `RAW_UQW1b_25`（仅去 DSA）：3 模均值 $\Delta E_1 \ge +0.015$ 且 `3/3 > 0`；
  3. 相对 `DSA_UQW1_25`（仅去 $B=300$）：3 模均值 $\Delta E_1 \ge +0.010$ 且 `3/3 > 0`；
  4. 相对 `RAW_UQW1b_100`（旧 $f_{\text{ref}}=100\text{ Hz}$）：3 模均值 $\Delta E_1 \ge +0.150$ 且 `3/3 > 0`；
  5. 过补偿 $m=0.5$ 上相对 `max(DSA_IQ70_25, DSA_W1b_25)`：`3/3 > 0` 且均值 $\ge +0.030$。
- **`G3_Synth_VolveTwin_U4_FullConfluence`**：在 Volve 全链路孪生体 `U4`（5 井 `[4–18 cyc/km]`）上：
  1. 未锚定 `RAW_UQW1_25` 均分低于 `RAW_R0`（复现实测掉分）；
  2. 统一 **`DSA_UQW1b_PAnc_25`** 均值与 $m=1.0$ 均超 `RAW_UQW1_25` 至少 `+0.050`，且均 $\ge \text{RAW\_R0}$；
  3. 过补偿 $m=0.5$ 上超过 `DSA_W1b_PAnc_25` 与 `DSA_IQ70_PAnc_25`。
- **`G4_Volve25_Confluence_Closure`**：在实测 Volve 4–25 Hz（`fld`）上：
  1. `RAW_R0` 与 Stage 4 存档偏差 $< 10^{-4}$；
  2. 上覆 `[1600, 2400] m` 画弧能量占比 `Smile_Frac` 较 `RAW_R0` 下降 $\ge 75\%$；
  3. 储层 `[2600, 3200] m` 谱斜率 $|\beta_{\text{res}}|$ 较 `DSA_R0_Eq` 与 `RAW_UQW1b_Unanc_25` 均下降 $\ge 50\%$；
  4. 主推 **`DSA_UQW1b_PAnc_25`**（及 `Amp` 伴随）在 `F5_R0band` 上超过 Stage 4 归档 `UQW1(25)` 至少 `+0.100`，且在宽带 `F5_Broad` 上超过 `RAW_R0`。
- **`G5_Volve45_Smile_And_Spectral_Closure`**：在实测 Volve 4–45 Hz（`fld45`）上：
  1. `RAW_R0` 与 `RAW_UQW1b_Unanc_100` 复现 Stage 4H 存档（误差 $< 10^{-3}$）；
  2. `DSA_UQW1b_PAnc_25` 将上覆 `[1600, 2400] m` `Smile_Frac` 较 Stage 4H `UQW1b`（`0.830`）降低 $\ge 70\%$；
  3. `DSA_UQW1b_PAnc_25 (m=1.0)` 的储层谱斜率 $|\beta_{\text{res}}|$ 较未锚定 `RAW_UQW1b_Unanc_100` 与 `DSA_R0_Eq` 均下降 $\ge 50\%$。
- **`G6_Volve45_Broadband_WellTie_And_Robustness`**：在实测 Volve 4–45 Hz `tie_fix`（8–40 Hz）上：
  1. **`DSA_UQW1b_PAnc_25`** 在全部 5 个 $m$ 上 $\ge \text{RAW\_R0}$ (`0.4740`)，且 $m=1.0$ 严格高于 `RAW_R0` 至少 `+0.030`；
  2. 在过补偿端 $m=0.6$ 严格超过 `DSA_W1b_PAnc_25` 与 `DSA_IQ70_PAnc_25`，且无单调崩盘（`tie(m=0.6) >= tie(m=1.4) - 0.030`）。
- **`G7_VelocityError_And_HighSNR_Safety`**：在 `U5`（高信噪比 15 dB + `-3%` 速度误差 + 横向 Q 梯度）上，`DSA_UQW1b_25` 相对 `RAW_UQW1_25` 满足 $\Delta E_1 \ge -0.015$ 且 $\Delta\text{SNR}_{\text{img}} \ge +0.5\text{ dB}$。

---

## 4. 诚实汇报与归档纪律
- 无论 `G1–G7` 全过还是局部未过，严格按 `judge_unified.py` 自动裁决原样落盘并写入 `RESULTS_ZH.md`，绝不事后修改代码、频段或门槛；特别地，对实测 `f_ref=25 Hz`、`f_ref=f (Amp)` 与 `f_ref=100 Hz` 的相位/深标差异（含零时移相关与最优深度滞后）做透明物理剖析。
- 实验完成后更新 `地震衰减补偿研究/README.md` 与 `目录标签.csv`，制作论文级高分辨率多联对比图，并按 `AGENTS.md` 提交推送至分支 `agy/uq-rtm-unified-20261124`。
