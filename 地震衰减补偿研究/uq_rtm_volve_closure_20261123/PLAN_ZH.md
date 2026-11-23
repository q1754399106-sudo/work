# 预注册实验方案：Volve 北海实测海上 OBC 地震资料闭环攻关（`uq_rtm_volve_closure_20261123`）

- **日期**：2026-11-23
- **工作区与分支**：`D:\ChatGPT\研究\.claude\worktrees\bosr-wsuq-20261123`（分支 `agy/bosr-wsuq-20261123`）
- **输入数据（严格只读）**：`D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\cache\fld`（124 炮、28,288 道、85 频 `4.00:0.25:25.00 Hz`、VTI-2008 迁移速度、名义 $Q=120$、缓存 `tfields.mat`、`ts.mat`、`state1.mat`；绝不修改原目录，绝不触碰 `fld45`）
- **测井盲测铁律**：算法全部参数（浅层锚定窗 `[1000, 1600] m`、参数化平滑频谱斜率 $\beta_{\text{anc}}$、盲结构张量倾角 $\phi_{\text{eff}}(x,z)$、Nyquist 防假频孔径上限 $\theta_{\max}=50^\circ$、沿用 Stage 4 既有的 $\pm 5\text{ Hz}$ 带内求和反极性噪声估计 `inv.mat`）**100% 仅由地震数据、采集几何与速度/$Q$ 模型决定，严禁在任何算法环节触碰 `well_F-5` 或 `well_F-14` 测井曲线**。

---

## 1. 实测物理诊断（纯地震观测，零测井先验）

对 Stage 4 实测缓存的纯地震物理审计揭示了原版 `fld_run.m`（及 `fld_post.m`）在实测 Volve 上全线落后于 `R0` 的三重叠加根因：

1. **相位参考频率错配（Stage 4 事后已识别）**：2008 迁移速度场属于地震频带内速度，`f_ref = 100 Hz` 引入系统性体时移；采用纯振幅补偿（`f_ref = f`）或带内中频参考（`f_ref = 14.5 Hz`）可消除体时移。
2. **未锚定的地表 $t^*$ 补偿叠加在强浅层蓝倾之上**：
   - `fld_prep.py` 直接将三维海上实测 OBC 记录（含 6 m 气枪震源鬼波 $\propto f^{+0.9}$ 与三维球面扩散相对二维的缺失半积分 $\propto f^{+0.5}$）喂入二维 Helmholtz 引擎并乘以 $(2\pi f)^2$。
   - 纯地震频谱审计测得：浅层 `[1000, 1600] m` 的 `R0` 频谱在 `[6, 22] Hz` 具有高达 **`+0.0541 / Hz`**（等效 $t^* = -17.2\text{ ms}$）的强蓝倾；即使深部储层 `[2600, 3200] m` 经历了真实地层衰减（浅至深实测 $\Delta t^*_{\text{emp}} \approx 7.7\text{–}10.9\text{ ms}$，$Q_{\text{eff}} \approx 77\text{–}108$），其 `R0` 频谱仍残余 **`+0.0298 / Hz`**（$-9.5\text{ ms}$）蓝倾。
   - 原版对地表全路径 $t^* \approx 23.3\text{ ms}$ 做补偿，将储层推向极端过蓝倾；而探索性逐频除以浅层经验谱 $S_{\text{anc}}(f)$ 又会把浅层 `6–14 Hz` 薄层干涉凹口反向注入储层，造成波形偏肥。
   - **本文方案——参数化平滑浅层斜率锚定（PSAE）**：在浅层 `[1000, 1600] m`、`x ∈ [-2500, 1500] m` 窗内对盲 `DSA_R0(f)` 振幅谱做 `[6, 22] Hz` 单参数对数线性拟合 $\ln S_{\text{anc}}(f) \approx a_0 + \beta_{\text{anc}} f$，构造无凹口平滑均衡因子：
     $$W_{\text{eq}}(f) = \exp\!\big(-\beta_{\text{anc}} (f - f_{\text{mid}})\big), \quad f_{\text{mid}} = 14.5\text{ Hz}$$
     并将走时衰减严格锚定至浅层基准：
     $$\Delta t^*_{\text{rel}}(x,z,s) = \max\!\big(t^*_{\text{dsa}}(x,z,s) - t^*_{\text{anc}},\, 0\big), \quad t^*_{\text{anc}} = \operatorname{median}_{[1000,1600]\text{ m}}\big(t^*_{\text{dsa}}\big)$$
3. **静止 OBC 全孔径宽角假频与 36.5% 出阵炮画弧**：
   - Volve 测线采用 50 m 炮间距、±5 km 超大偏移距静止海底电缆（OBC），上覆地层 $v \approx 2000\text{ m/s}$ 在 $f_{\max}=25\text{ Hz}$ 下的炮采样防假频临界角为 $\theta_{\text{alias}} = \arcsin\!\left(\frac{2000}{2\times 50\times 25}\right) \approx 53^\circ$，且储层窗内 36.5% 炮-点对的镜面反射点落在检波阵外，导致 `1600–2400 m` 出现密集倾斜画弧假象，且原版按道反极性估噪把相干画弧误判为信号（`SNR = 12–20 dB` 虚高）。
   - **本文方案——叠前盲倾角镜面孔径 + Nyquist 防假频角限（`DSA`, $\theta_{\max} = 50^\circ$）**：
     - 盲结构张量倾角 $\phi_{\text{eff}}(x,z)$ 从原始未补偿 `R0` 估计（高斯平滑 $\sigma = 60\text{ m}$，限幅 $\pm 30^\circ$，乘以相干度 $c_{\text{st}}$）；
     - 镜面孔径权重 $w_{\text{dsa}}(x,z,s)$ 同时施加阵内余弦过渡（$d_{\text{in}} \in [-150, 250]\text{ m}$）与源/检角度防假频余弦过渡（$|\theta_s|, |\theta_r| \in [40^\circ, 50^\circ]$，下限 $10^{-4}$）；
     - $t^*_{\text{dsa}}(x,z,s)$ 直接利用缓存的每炮单程 eikonal $t^*$ 场（`tfields.mat`）按倾角修正镜面点 $x_{\text{mir}}$ 插值构造；
     - 噪声功率比 `iv(x,z,f)` 严格沿用 Stage 4 已缓存的 $\pm 5\text{ Hz}$ 带内求和反极性估计 `cache/fld/inv.mat`（3 频冒烟测试发现单频分炮翻转会破坏非镜面炮驻相抵消导致高频假性 `-9 dB` 窒息补偿，故锁定 Stage 4 同口径带内 `inv.mat`，确保与 Stage 4 对比完全苹果对苹果，所有改进纯粹来自 `DSA` 与 `PSAE`）。

---

## 2. 单遍 85 频全波叠前 2×2 因子对照设计

在同一遍全 85 频（`4.00:0.25:25.00 Hz`）Helmholtz 求解循环中，保存每频实部成像切片后统一应用盲估 $W_{\text{eq}}(f)$，完整产出 9 组对照（$m \in \{0.5, 0.7, 1.0, 1.4, 2.0\}$，名义 $Q=120$）：

| 编号 | 方法代号 | 孔径屏蔽 | 频谱锚定与 $t^*$ | 相位设置 | 角色 |
|---|---|---|---|---|---|
| 1 | `RAW_R0` | 全炮 (`RAW`) | 无 | — | 精确复现 Stage 4 原始 `R0`（校验误差 $<10^{-4}$） |
| 2 | `DSA_R0` | `DSA` ($\theta_{\max}=50^\circ$) | 无 | — | 剥离纯 `DSA` 去假频/去出阵炮效应 |
| 3 | `DSA_R0_Eq` | `DSA` | `PSAE` ($W_{\text{eq}}$)，无 $Q$ 补偿 | — | **物理公平的浅层均衡未补偿基线**（浅平、深衰） |
| 4 | `RAW_UQW1_Unanc` | `RAW` | 无（全地表 $t^*$） | 纯振幅 (`fref=f`) | 复现 Stage 4 `UQW1_amp`（复用同口径 `inv.mat`） |
| 5 | `DSA_UQW1_Unanc` | `DSA` | 无（全地表 $t^*$） | 纯振幅 (`fref=f`) | 剥离“只加 `DSA` 不做浅层锚定” |
| 6 | `RAW_UQW1_PAnc` | `RAW` | `PSAE` + $\Delta t^*_{\text{rel}}$ | 纯振幅 (`fref=f`) | 剥离“只做浅层锚定不加 `DSA`” |
| 7 | `DSA_IQ70_PAnc` | `DSA` | `PSAE` + $\Delta t^*_{\text{rel}}$ | 纯振幅 (`fref=f`) | 公平竞争者：反 Q 70 dB |
| 8 | `DSA_W1_PAnc` | `DSA` | `PSAE` + $\Delta t^*_{\text{rel}}$ | 纯振幅 (`fref=f`) | 公平竞争者：单点 Wiener (`W1`) |
| 9 | **`DSA_UQW1_PAnc`** | **`DSA`** | **`PSAE` + $\Delta t^*_{\text{rel}}$** | **纯振幅 (`fref=f`)** | **核心主方法（Hero）** |
| 10 | `DSA_UQW1_PAnc_Disp15` | `DSA` | `PSAE` + $\Delta t^*_{\text{rel}}$ | 带内频散 (`fref=14.5 Hz`) | 因果频散对照变体 |

---

## 3. 预注册评价指标与硬性通过门限（`judge_volve_closure.py`）

1. **井震标定口径（与原 Stage 4 `fld_eval.m` 逐位对齐 + 宽带扩展）**：
   - **`F5_R0band`**：严格复现 `fld_eval.m`：`lo = 3.515625 cyc/km`, `hi = 12.890625 cyc/km`，40 m 垂直平滑高通，井旁 ±2 道（5 道均值），允许滞后 $L \in [-3, 3]$ 采样点（$\pm 30\text{ m}$）。原 Stage 4 基线：`R0 = 0.6259`，`UQW1_amp = [0.5032, 0.5289, 0.5532, 0.5723, 0.5882]`（均值 `0.5492`）。
   - **`F5_Broad`**：固定宽带 `[3.515625, 16.0] cyc/km`（检验超越 `R0` `-10 dB` 上界的高波数真实反射率恢复）。
   - **`F14_R0band`**：离线 190–540 m 参考井，仅作描述记录，不入硬门限。
2. **纯地震物理指标**：
   - **`Res_Slope`**：储层窗 `[2600, 3200] m` 在 `[6, 22] Hz` 的对数线性频谱斜率 $\beta_{\text{res}}$（理想无偏平谱为 `0.000 / Hz`）。
   - **`Smile_Frac`**：上覆假频区 `[1600, 2400] m, x ∈ [-2500, 1500] m` 二维 $k_x\text{–}k_z$ 谱中陡倾假频/画弧能量占比（$|k_x / k_z| > \tan 25^\circ$）。
   - **`Res_Cen`**：储层窗深度波数谱质心（cyc/km）与 `-10 dB` 带宽。

### 5 项预注册自动判决门限
- **`V1_Beats_Stage4_UQW1`**：`DSA_UQW1_PAnc` 在 `F5_R0band` 上的 5-`m` 均值、`m=0.5`、`m=1.0` 全部超过原 Stage 4 归档 `UQW1_amp` 至少 **`+0.050`**。
- **`V2_Beats_Uncomp_Equalized`**：在浅层均衡公平基线上，`DSA_UQW1_PAnc` 在 `m=1.0` 和 5-`m` 均值上均优于 `DSA_R0_Eq`（`F5_R0band` 与 `F5_Broad` 两者中至少 `F5_Broad` 严格提升且 `F5_R0band` 均值与 `m=1.0` 均不低于 `DSA_R0_Eq`）。
- **`V3_Beats_Original_RAW_R0`（实测终极闭环门槛）**：`DSA_UQW1_PAnc` 在 `m=1.0` 或 5-`m` 均值上达到或超过原 Stage 4 `RAW_R0`（`F5_R0band >= RAW_R0` 且 `F5_Broad >= RAW_R0`）。
- **`V4_Overcomp_Robustness_VSP`**：在 VSP 低 Q 过补端 `m=0.5`，`DSA_UQW1_PAnc` 在 `F5_R0band` 上超过 `DSA_IQ70_PAnc` 和 `DSA_W1_PAnc`，且 `tie(m=0.5) >= tie(m=2.0) - 0.03`（无单调过补崩盘）。
- **`V5_Seismic_Physics_Closure`**：(a) `DSA_UQW1_PAnc (m=1.0)` 的储层频谱斜率绝对值 $|\beta_{\text{res}}|$ 较 `DSA_R0_Eq` 与 `RAW_UQW1_Unanc` 均下降 $\ge 50\%$；(b) `DSA_R0` 较 `RAW_R0` 降低上覆 `[1600, 2400] m` 陡倾画弧能量占比 `Smile_Frac`。

---

## 4. 失败处理与诚实汇报原则
- 无论 `V1–V5` 全部通过还是部分未过，均按 `judge_volve_closure.py` 自动裁决原样落盘 `results/volve_closure_verdict.json` 并写入 `RESULTS_ZH.md`，绝不事后修改脚本、频段或门限。
- 全程严禁将任何 `.mat` / 频率切片缓存提交到 Git 仓库（全部存放于工作区外 `scratch/`）。
