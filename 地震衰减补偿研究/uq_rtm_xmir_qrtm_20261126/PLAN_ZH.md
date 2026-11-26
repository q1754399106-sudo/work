# 预登记实验方案：`uq_rtm_xmir_qrtm_20261126`
## 升级弯曲射线 `xmir_bent` 与引入黏声波方程对照、硬刚真 Q-RTM 全量盲测方案

- **创建日期**：2026-11-27
- **实验分支**：`agy/uq-rtm-xmir-qrtm-20261126`（自 `origin/main` 切出）
- **科研诚信约束**：严格执行全局六条铁律与 Gate `X0` 机械审计（`scripts/audit_report_integrity.py` 退出码 `0` 零违规）；全部门控门槛锁定在 `scripts/judge_xmir_qrtm.py`，开跑前记录 SHA-256 至 `PRE_RUN_HASHES.txt`，开跑后禁止篡改。

---

## 1. 科学动机与两大核心升级

1. **升级 `xmir`（从均匀直射线近似升级为二维非均匀弯曲射线 `xmir_bent`）**：
   - 既有 `DSA` 使用均匀介质直射线公式 $x_{\text{mir,str}} = x + \Delta z_r \tan(\theta_s^{\text{chord}} + 2\phi)$，在强速度梯度（如 Volve 海底 `1687 m/s` 至储层 `3956 m/s`）与倾角地层中忽略 Snell 射线弯曲，导致局部入射角与镜面出射点偏差达数百米至上千米。
   - **升级方案 (`build_dsa_bent`)**：
     1. 由已求解的程函走时场梯度 $\nabla T_s(\mathbf{x}, k)$ 直接提取精确弯曲射线局部入射角 $\theta_s^{\text{eik}} = \operatorname{atan2}(\partial_x T_s, \partial_z T_s)$；
     2. 由局部倾角 $\phi(\mathbf{x})$ 给出精确 Snell 反射起飞角 $\theta_r = \theta_s^{\text{eik}} + 2\phi$（并检测背向入射与下折转弯射线）；
     3. 在完整二维非均匀慢度场 $s(z, x) = 1/v_{\text{mig}}(z, x)$ 中执行逐层同步二阶 Runge–Kutta（RK2）上行运动学射线追踪（大网格采用炮无关 3D RK2 查找表），获得精确弯曲射线出射点 $x_{\text{mir,bent}}$、地表出射角 $\theta_{r,\text{surf}}$、转弯/受困射线掩码及倾角自洽双程衰减 $t^*_{\text{bent}}$。
2. **引入黏声波方程真 Q-RTM 对照体系（硬刚真 Q-RTM，拒绝稻草人基线）**：
   - 在完全一致的 Kjartansson/Kolsky 因果复速度频率域黏声波方程（`helm_op_visco`：物理内部严格复速度、C-PML 边界严格实数无损吸收）中实现 **8 种真黏声波方程偏移对照**（外加 1 种波方程+贝叶斯融合探索）：
     - `WE_Q_Adj`（标准伴随黏声 RTM：精确频散、双倍衰减）
     - `WE_Q_DispOnly` / `DSA_WE_DispOnly`（真黏声纯频散波方程 RTM：RAW 与 +Bent-DSA）
     - `WE_QRTM_Unreg`（无正则化真黏声 Q-RTM：无噪真 $Q$ 理论金标）
     - `WE_QRTM_ZhuLP`（Zhu et al. 2014 频率余弦截断稳定真 Q-RTM：全频频散 + 高频衰减补偿平滑截断）
     - `WE_QRTM_Reg` / `DSA_WE_QRTM_Reg`（Wang/Sun–Zhu 2018 增益受限 Tikhonov 稳定真 Q-RTM：RAW 与 +Bent-DSA）
     - `DSA_WE_QRTM_W1b`（**最强确定性真 Q-RTM 对手**：真黏声波方程 + 同款 `xmir_bent` 门控 + 同款 `B=300` 噪声自适应！）
     - `WE_Disp_UQW1b_25`（融合探索：真黏声频散波方程 + 有限频波方程 $t^*_{\text{WE}}$ + 贝叶斯 `UQW1b` 振幅核）

---

## 2. 四大工作包（WP1–WP4）与全新盲测矩阵

| 工作包 | 数据集 / 模型 | 种子与信噪比校准 | 核心检验问题 |
|---|---|---|---|
| **WP1** | `wp1_physics_verify` | 解析 $v(z)=1500+0.8z$ + 无噪气烟囱 + Volve 10m 速度场 | 验证 `xmir_bent` 相对解析圆弧 Snell 真解降噪 $\ge 80\%$；验证真 Q-RTM 无噪恢复且 C-PML 稳定；审计 Volve 实测直/弯射线几何差异 |
| **WP2** | `S1_Chimney_Base`、`S2_StrongGrad_SteepDip`、`S3_SharpQ_Multipath`、`S4_NonScalar_SpatialQ` | 新种子 `20261501–20261504`；开跑前冒烟校准至 `RAW_R0 SNR_img ∈ [+2, +6] dB`；`l30` 50m 相关噪 | 全 16 法在无噪 $m=1$ 与含噪多 $m$ 下正面硬刚：基线气烟囱（`S1`）、强梯度陡倾非对称缆（`S2`）、尖锐低 $Q$/低速多路径透镜（`S3`，主动寻找真 Q-RTM 胜出边界）、非标量错位+深浅反向 $Q$ 误差（`S4`） |
| **WP3** | `S5_FD_GSLS_NoInvCrime` | 新种子 `20261505`；独立 **时间域 8 阶交错网格 GSLS 有限差分 (`fd2d_visco.m`, 4 弛豫机制)** 正演 | 彻底消除“反演犯罪”质疑：在完全独立时间域 GSLS 正演数据上复验全 16 法硬刚真 Q-RTM |
| **WP4** | `volve25_xmir_qrtm` | Volve 北海 OBC 4–25 Hz（**全 85 频 $\Delta f = 0.25\text{ Hz}$ 零混叠**），统一 $m \in \{0.5, 1.0, 2.0\}$ | 全 20 法实测硬刚：直射线 `xmir_str` vs 弯曲射线 `xmir_bent` vs 7 种真黏声波方程 Q-RTM（含原始无锚定与公平 `PSAE` 浅层锚定版本） |

---

## 3. 预登记机械门控（冻结于 `scripts/judge_xmir_qrtm.py`，不得事后修改）

- **`X0_Integrity_Gate`**：`audit_report_integrity.py` 退出码 `0`（全行、全列、反超核查、高频活跃增益 `G_HF`、反夸大词、代码逐行接线核查全部零报错）。
- **`X1_Analytical_BentRay_xmir`**：解析圆弧测试中 4 个非零倾角（`±10°, ±20°`）`xmir_bent` 中位误差降幅均 $\ge 80\%$，且落盘 Volve 10m 实测几何偏差审计。
- **`X2_TrueQRTM_NoStrawman_Proof`**：`S1–S3` 无噪 `m=1` 下 `WE_QRTM_Unreg` 在 3/3 模型上均实现 `E1 >= RAW_R0 + 0.04` 且 `E_amp` 较 `RAW_R0` 压降 $\ge 60\%$（证明真 Q-RTM 对手非稻草人）。
- **`X3_BentRay_Imaging_Benefit_S2`**：强梯度+陡倾+非对称缆 `S2` 上（`RAW_R0 SNR_img ∈ [+2, +6] dB`），`DSA_Bent_UQW1b_25` 较 `DSA_Str_UQW1b_25` 实现 `E1_mean >= +0.005` 或 `E_amp` 降低 $\ge 10\%$，且 `Smile_Frac <= 12%`。
- **`X4_HardHeadToHead_UncertainQ`**：`S1–S4` 全 4 模型 `RAW_R0 SNR_img ∈ [+2, +6] dB`，在过补偿高压端（`S1–S3` $m=0.5$、`S4` $m=0.6$），`DSA_Bent_UQW1b_25` 在 **4/4** 模型上高于全部 3 种标准确定性真 Q-RTM（`WE_QRTM_Unreg`, `WE_QRTM_ZhuLP`, `WE_QRTM_Reg`），且在 **$\ge 3/4$** 模型上高于或持平（$\ge -0.005$）最强噪声自适应真 Q-RTM（`DSA_WE_QRTM_W1b`）。
- **`X5_Multipath_Crossover_Audit_S3`**：尖锐多路径透镜 `S3` 上如实定量刻画无噪 $m=1$（真 Q-RTM 物理优势区）与含噪不确定 $Q$（贝叶斯 UQ 稳健区）的精确交叉边界，且 `DSA_Bent_UQW1b_25` 与 `WE_Disp_UQW1b_25` 最坏 `min E1 >= RAW_R0 - 0.010`。
- **`X6_Independent_GSLS_FD_S5`**：独立时间域 8 阶 GSLS 有限差分 `S5`（`RAW_R0 SNR_img ∈ [+2, +6] dB`）上，`DSA_Bent_UQW1b_25` 在 5-$m$ `E1_mean` 上高于 `RAW_R0`、`WE_QRTM_Unreg`、`WE_QRTM_Reg`，在 $m=0.5$ 过补偿端高于 `DSA_WE_QRTM_W1b`，且 `E_amp` 较 `WE_Q_DispOnly` 压降 $\ge 25\%$。
- **`X7_Volve25_Field_Showdown`**：Volve 4–25 Hz 实测全 85 频零混叠下，`DSA_Bent_UQW1b_PAnc_25` 真实开机（`G_HF[18–24Hz] >= 1.20x` vs `DSA_Bent_R0_Eq`）、`Smile_Frac <= 8%`、过补偿 $m=0.5$ 稳守 `F5_R0band >= RAW_R0 - 0.010`，并 100% 全量呈现与 `DSA_Str` 及 7 种真黏声 Q-RTM 的实测对比与算力账本。
