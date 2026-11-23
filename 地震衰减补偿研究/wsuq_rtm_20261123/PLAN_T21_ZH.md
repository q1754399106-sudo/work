# T2.1 预注册验证方案 (`wsuq_rtm_20261123`)：真实程函多炮基准与全波动方程 2D 黏声气烟囱 RTM 对比

- **前置门控**：T2.0 (`85c2d5a`, `RESULTS_STEP0_ZH.md`) 已通过 **G0-1**（全部 3 条边界恒等式误差 $\le 10^{-6}$）与 **G0-2**（Volve F-5 $R_{\text{CS}}(35\text{ Hz}, m=1) = 1.35$, $R_{\text{CS}}(22\text{ Hz}, m=0.5) = 1.60$, $R_{\text{CS}}(45\text{ Hz}, m=0.5) = 5.76 / +7.6\text{ dB}$）。
- **目的**：在严格控制宿主机内存安全（峰值 $< 0.8\text{ GB}$，与同机运行的其他 MATLAB 任务零冲突）的前提下，以两条正交路径验证 `WS-UQ1` / `WS-VA-UQ1` 相对 `UQW1_old` / `VA_old` / `IQ70_old` / `R0` 的非零和性能增益：
  1. **Part A（`scripts/wsuq_eikonal_bench.m`，零 LU，约 30 s）**：直接读取真实 Volve F-5（124 炮 `TS_abs` / `TS_anc227` + 真实 `well_F-5.bin` 反射系数）与 `D2`（51 炮 `TS_abs` / `TS_anc` + 真实弹性参考 `ref.R0`），在全宽频 `4–45 Hz`、双噪声水平（10 dB / 3 dB）与 5 档 $Q$ 失配 $m \in \{0.5, 0.7, 1.0, 1.4, 2.0\}$ 下对比全部算子。
  2. **Part B（`scripts/wsuq_rtm2d_bench.m`，全 2D 黏声 Helmholtz 波动方程正演 + 2D 频域多炮 RTM，约 4 min，峰值内存 $< 0.6\text{ GB}$）**：
     - 零反演犯罪：正演使用精确复速度 $v_c(z,x,f)$ 黏声 Helmholtz（含 4 层沉积介质 + 强低 $Q$ 气烟囱 $Q_{\text{chim}}=22$ + 3 个深部背斜/平层目标反射面，内部网格 `140 x 240`, `dx = 8 m`, `27 炮`, `f = 4:0.5:40 Hz` 共 73 频点）；
     - 偏移仅使用平滑背景速度 $v_{\text{band}}(25\text{ Hz})$ + FMM 程函双程 $t_s^*$ 与盲图像域极性翻转噪声估计 `iv`；
     - 在低噪（`clean`, 15 dB）与现实强噪（`noisy`, 3 dB）两种数据信噪比下，使用 **4 项非饱和指标**（宽频相关 `E1`、高频波数相关 `E1_HF`、烟囱内外振幅恢复误差 `E_amp`、图像域信噪比 `SNR_img_dB`）全面对比 `R0`, `IQ70_old`, `UQW1_old`, `VA_old`, `WS_UQ1`, `WS_VA_UQ1`。

## 预注册机械判定准则（`scripts/judge_t21.py`）

- **`T21_Eik_Volve`**：在 Volve F-5 真实 124 炮 `TS_abs`（`4–45 Hz`, 3 dB 档）上，`WS_UQ1` 跨 $m$ 平均输出 `SNR_dB` 优于 `UQW1_old` $\ge +1.0\text{ dB}$，且 $m=0.5$ 处优势 $\ge +2.5\text{ dB}$，高频相关 `E1_HF` 同步提升。
- **`T21_RTM2D_Gain`**：在真实 2D 全波气烟囱 RTM（`noisy` 档）烟囱下窗 `W_chim` 上，`WS_UQ1` 相对 `UQW1_old` 满足：跨 $m$ 平均 `E1_HF` 严格更高、平均 `SNR_img_dB` 提升 $\ge +0.5\text{ dB}$，且最恶劣过补偿点 $m=0.5$ 的 `E1` 严格更高。
- **`T21_RTM2D_NoHarm`**：在 `clean` 档准 $Q$ ($m=1.0$) 上，`WS_UQ1` 的 `E1` 不低于 `UQW1_old` 超过 `0.01`。
