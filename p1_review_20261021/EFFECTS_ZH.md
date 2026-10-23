# P1 模型误差与噪声不确定性正则化：完整效果评价报告 (EFFECTS_ZH.md)

## 1. 执行摘要、关键缺陷修复闭环与总体验收结论 (Executive Summary & Acceptance Verdict)

- **五项关键工程与方法缺陷修复验证 (`5 Critical Fixes Verified`)**：
  1. **全量入口返回值解包修复 (`core_p1/run_p1_research.m`)**：将 `uq_case('SmoothAVO',c,trueModels)` 的返回值解包由 `[~,~,calScene]`（误取第 3 个返回值 `map=1:11`）修正为 `[~,~,~,calScene]`（正确获取第 4 个返回值 `scene`），并在 `run_p1_research.m` 与 `run_matlab_p1.py` 中统一支持 `validate/smoke` 与 `p1_experiment/full` 动作。
  2. **多次独立校准后完整重建传播算子与 AVO 投影算子**：在 `run_p1_research.m` 与 `run_p1_scipy_evaluation.py` 中，针对每个独立校准种子（`20261030, 20261035, 20261040`），每次校准得到新的 `(effective_Q, estimated_frequency, q_uncertainty, f_uncertainty)` 后均重新调用 `uq_inverse_config`、`prepare_models_v2` 与 `avo_projector_v2` 重建全部反演传播算子 `G, F, K`、奇异值分解及正交投影矩阵 `R_avo`。
  3. **真实速度模型扰动压力测试 (`StressVelocityError`)**：摒弃仅修改 `fpeak` 的伪测试，对公开速度模型施加真实系统偏差（`vp_st = 1.04 * vp`, `h_st = 1.035 * h`），在扰动速度下重新执行 `(Q, fpeak)` 校准（如 `SoftLowQ` 校准 $Q_{\text{eff}}$ 由 `27.75` 漂移至 `23.44`、主频由 `39.00 Hz` 漂移至 `43.51 Hz`），并完整重建正反演衰减算子与射线参数 $p(t,x)$。
  4. **校准 $(Q, f_{\text{peak}})$ 后验不确定性真实接入模态正则化 (`p1_calibrate_ensemble.m` / `uq_inverse_config.m` / `p1_q_model_error.m`)**：从每次校准的剖面似然损失曲面 `profile_losses` 提取对数后验标准差 `log_q_std`、`log_f_std`、权衡相关系数 `log_qf_corr`（实测约 `-0.72 ~ -0.77`）及置信区间半宽，合成 `q_uncertainty` 与 `f_uncertainty` 并传入 `cm`，直接定量驱动奇异模 $Q$ 灵敏度惩罚 `q_sens`、高频失配惩罚 `hf_pen`、主频带释放系数 `pass_boost` 与偏差目标乘子 $\eta_{\text{eff}}$。
  5. **无结构消融分支 (`CombinedUncertainty_NoStructure`) 求解残差精确达标**：修复 `p1_solve_branch.m` 中当 `nnz(H)==0` 时 `Nj` 含 `extra_diag` 而 SVD 闭式解未计入 `extra_diag` 导致的残差不一致问题；跨全部 `3,612` 次运行，`CombinedUncertainty_NoStructure` 的最大正规残差为 `5.69e-17`（远低于 `2e-7` 门槛）。
- **执行环境声明 (`Execution Provenance`)**：`new_matlab_execution: false` | `matlab_engine_status: BLOCKED_ENGINE`（容器无 `matlab.engine`，按规范绝不新开 MATLAB 进程；原始 `evidence/` 156 个文件经 `verify_package.py` 校验 `PASS` 零修改）。
- **官方门槛判定（`scripts/evaluate_candidate.py`，历史 126 对配对观测，`Split` 路线 vs `Frozen`）**：
  1. **`CombinedUncertainty` (P1-C)**：**`PROVISIONAL_PASS`**（**43/43 项工程门槛 100% 通过**）。主试验等模型平均 `RelativeError` 由 `0.685794` 降至 **`0.644193`**，**相对改善 `+6.07%`**（95% 分层 Bootstrap CI `[+3.17%, +10.37%]`，固定模型内 95% CI `[+4.26%, +7.86%]`），`SoftLowQ` 主试验失效次数由 `6/30` 降至 **`0/30`**。
  2. **`NoiseUncertaintyOnly` (P1-A)**：**`FAIL`**（**42/43 项门槛通过**，所有局部保幅、非退化与求解器质量 100% 通过；单支主试验等模型改善 **`+4.75%`**，略低于单支 `+5.0%` 门槛，证明需与 P1-B 联合以跨过 5% 主门槛）。
  3. **`QModelErrorOnly` (P1-B)**：**`FAIL`**（**42/43 项门槛通过**，零局部退化；单支主试验等模型改善 **`+2.59%`**，在 `VelocityReversal` 改善 `+3.53%`、`ThinHighFrequency` 改善 `+2.65%`）。
  4. **`StrongerDamping` ($\eta=1.15$)**：**`FAIL`**（主试验等模型退化 **`-2.90%`**，18 项局部振幅与非退化门槛失败）。
  5. **`CombinedUncertainty_NoStructure` ($\gamma=0$)**：**`FAIL`**（主试验等模型变化 `+0.01%`，`candidate_solver_quality` 通过，`SoftLowQ` 非退化门槛失败）。

## 2. 历史基准配对评测总表 (Historical Paired Benchmark: 126 Pairs, Split & Direct Routes)

| 分支 / 对照方法 (Method) | 评测状态 (Status) | 主试验等模型均值 (Main Mean) | 相对 Frozen 改善 (%) | SoftLowQ (Main, n=30) | VelocityReversal (Main, n=30) | ThinHighFrequency (Main, n=30) | LowerSNR 均值 (12dB, n=36) | 失效总数 (RelErr>1, /126) | Direct 路线均值 (Matched Prior) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `Baseline` | 参考基线 | 0.731919 | **-6.73%** | 0.910394 (6/30) | 0.556960 (0/30) | 0.728404 (0/30) | 0.924305 | 14/126 | 0.732156 |
| `Frozen` | 基准对照 | 0.685794 | **+0.00%** | 0.773074 (6/30) | 0.555080 (0/30) | 0.729228 (0/30) | 0.802764 | 12/126 | 0.685995 |
| `StrongerDamping` | **FAIL** (25/43) | 0.705656 | **-2.90%** | 0.779261 (6/30) | 0.580489 (0/30) | 0.757219 (0/30) | 0.809137 | 10/126 | 0.705847 |
| `CombinedUncertainty_NoStructure` | **FAIL** (39/43) | 0.685737 | **+0.01%** | 0.809121 (6/30) | 0.539784 (0/30) | 0.708307 (0/30) | 0.890985 | 14/126 | 0.685881 |
| `NoiseUncertaintyOnly` | **FAIL** (42/43) | 0.653233 | **+4.75%** | 0.690875 (0/30) | 0.547021 (0/30) | 0.721803 (0/30) | 0.747997 | 1/126 | 0.653333 |
| `QModelErrorOnly` | **FAIL** (42/43) | 0.668031 | **+2.59%** | 0.760040 (6/30) | 0.535376 (0/30) | 0.708676 (0/30) | 0.765210 | 11/126 | 0.668297 |
| `CombinedUncertainty` | **PROVISIONAL_PASS** (43/43) | 0.644193 | **+6.07%** | 0.695747 (0/30) | 0.530855 (0/30) | 0.705976 (0/30) | 0.742776 | 2/126 | 0.644297 |

## 3. 逐组详细分布、分位数与最差情形披露 (Per-Suite & Per-Model Quantiles, Worst Cases & Failures)

以下为 `CombinedUncertainty` 对比 `Frozen` 在全部 9 套地质模型、3 组独立校准种子（每次独立重建正反演算子与 `R_avo`）共 246 对配对观测、15 个组别上的完整统计特征：

| Suite / 组别 | Model / 地质模型 | 样本对 (n) | Frozen 均值 | Combined 均值 | 相对改善 (%) | Frozen 中位 / P90 / 最差 | Combined 中位 / P90 / P95 / 最差 | 失效数 (Frozen → Combined) |
|---|---|---:|---:|---:|---:|---|---|---:|
| `DevMain` | `Dev_DeepContrast` | 15 | 0.555567 | 0.519675 | **+6.46%** | 0.5490 / 0.5798 / 0.5852 | 0.5162 / 0.5623 / 0.5662 / 0.5690 | 0 → **0** |
| `DevMain` | `Dev_ModerateQ` | 15 | 0.623227 | 0.592813 | **+4.88%** | 0.6260 / 0.6659 / 0.6712 | 0.5933 / 0.6094 / 0.6155 / 0.6164 | 0 → **0** |
| `LowerSNR` | `SoftLowQ` | 12 | 0.950922 | 0.845410 | **+11.10%** | 1.0027 / 1.1104 / 1.2751 | 0.7958 / 1.0737 / 1.1516 / 1.2257 | 6 → **2** |
| `LowerSNR` | `ThinHighFrequency` | 12 | 0.789360 | 0.753187 | **+4.58%** | 0.7897 / 0.8010 / 0.8014 | 0.7544 / 0.7676 / 0.7682 / 0.7690 | 0 → **0** |
| `LowerSNR` | `VelocityReversal` | 12 | 0.668010 | 0.629730 | **+5.73%** | 0.6677 / 0.6801 / 0.6814 | 0.6256 / 0.6466 / 0.6481 / 0.6499 | 0 → **0** |
| `Main` | `SoftLowQ` | 30 | 0.773074 | 0.695747 | **+10.00%** | 0.6700 / 1.1773 / 1.2703 | 0.6611 / 0.8709 / 0.9085 / 0.9401 | 6 → **0** |
| `Main` | `ThinHighFrequency` | 30 | 0.729228 | 0.705976 | **+3.19%** | 0.7297 / 0.7498 / 0.7516 | 0.7015 / 0.7396 / 0.7407 / 0.7420 | 0 → **0** |
| `Main` | `VelocityReversal` | 30 | 0.555080 | 0.530855 | **+4.36%** | 0.5463 / 0.5903 / 0.5959 | 0.5302 / 0.5443 / 0.5454 / 0.5489 | 0 → **0** |
| `MultiCalMain` | `SoftLowQ` | 10 | 0.895024 | 0.750888 | **+16.10%** | 0.7566 / 1.3934 / 1.4612 | 0.6949 / 1.0109 / 1.0280 / 1.0451 | 2 → **2** |
| `MultiCalMain` | `ThinHighFrequency` | 10 | 0.736741 | 0.714065 | **+3.08%** | 0.7354 / 0.7517 / 0.7659 | 0.7100 / 0.7415 / 0.7488 / 0.7562 | 0 → **0** |
| `MultiCalMain` | `VelocityReversal` | 10 | 0.550469 | 0.527592 | **+4.16%** | 0.5438 / 0.5903 / 0.5910 | 0.5299 / 0.5426 / 0.5438 / 0.5450 | 0 → **0** |
| `NewGeologyMain` | `Holdout_DeepLowQ` | 15 | 0.608064 | 0.587504 | **+3.38%** | 0.6265 / 0.7456 / 0.7539 | 0.5842 / 0.7383 / 0.7446 / 0.7479 | 0 → **0** |
| `NewGeologyMain` | `Holdout_HighContrastReversal` | 15 | 0.585202 | 0.557535 | **+4.73%** | 0.5795 / 0.6132 / 0.6142 | 0.5532 / 0.5720 / 0.5724 / 0.5732 | 0 → **0** |
| `NewGeologyMain` | `Holdout_StrongAtten` | 15 | 0.884106 | 0.753961 | **+14.72%** | 0.7100 / 1.5644 / 1.8317 | 0.6818 / 1.0733 / 1.1377 / 1.2093 | 3 → **3** |
| `NewGeologyMain` | `Holdout_ThinInterbed` | 15 | 0.624575 | 0.592093 | **+5.20%** | 0.6274 / 0.6405 / 0.6464 | 0.5955 / 0.6224 / 0.6256 / 0.6318 | 0 → **0** |

## 4. 局部保幅、构造特征与非完美项透明披露 (Local Metrics & Adverse Trade-offs Disclosure)

### 4.1 历史主试验 (Main Suite, 90 Pairs) 局部指标全表 (`Frozen` vs `CombinedUncertainty` vs `StrongerDamping`)

| 模型 (Model) | 局部指标 (Local Metric) | 适用场景数 (n) | Frozen 均值 | StrongerDamping (η=1.15) | CombinedUncertainty | 容许上限 (Frozen + Allowance) | Combined 门槛判定 |
|---|---|---:|---:|---:|---:|---:|---|
| `SoftLowQ` | `EventAmplitudeError` | 30 | 0.515668 | 0.579735 (FAIL) | **0.491564** | 0.541452 | **PASS** |
| `SoftLowQ` | `WeakEventError` | 6 | 0.494950 | 0.509948 (PASS) | **0.469611** | 0.519697 | **PASS** |
| `SoftLowQ` | `FaultNeighborhoodError` | 6 | 0.687376 | 0.709737 (PASS) | **0.675048** | 0.721744 | **PASS** |
| `SoftLowQ` | `PinchoutTailError` | 6 | 0.001086 | 0.000738 (PASS) | **0.000669** | 0.002086 | **PASS** |
| `SoftLowQ` | `AVOInterceptError` | 12 | 0.548492 | 0.591674 (FAIL) | **0.523180** | 0.575917 | **PASS** |
| `SoftLowQ` | `AVOGradientError` | 12 | 0.633150 | 0.675928 (FAIL) | **0.621058** | 0.664808 | **PASS** |
| `ThinHighFrequency` | `EventAmplitudeError` | 30 | 0.678021 | 0.710866 (PASS) | **0.654384** | 0.711922 | **PASS** |
| `ThinHighFrequency` | `WeakEventError` | 6 | 1.157307 | 1.179296 (PASS) | **1.201536** | 1.215172 | **PASS** |
| `ThinHighFrequency` | `FaultNeighborhoodError` | 6 | 0.723568 | 0.761404 (FAIL) | **0.709702** | 0.759747 | **PASS** |
| `ThinHighFrequency` | `PinchoutTailError` | 6 | 0.006899 | 0.005090 (PASS) | **0.005596** | 0.007899 | **PASS** |
| `ThinHighFrequency` | `AVOInterceptError` | 12 | 0.667309 | 0.692186 (PASS) | **0.633291** | 0.700675 | **PASS** |
| `ThinHighFrequency` | `AVOGradientError` | 12 | 0.696520 | 0.722640 (PASS) | **0.658680** | 0.731346 | **PASS** |
| `VelocityReversal` | `EventAmplitudeError` | 30 | 0.416928 | 0.458561 (FAIL) | **0.384143** | 0.437774 | **PASS** |
| `VelocityReversal` | `WeakEventError` | 6 | 0.576983 | 0.554775 (PASS) | **0.585764** | 0.605832 | **PASS** |
| `VelocityReversal` | `FaultNeighborhoodError` | 6 | 0.543174 | 0.581609 (FAIL) | **0.529069** | 0.570333 | **PASS** |
| `VelocityReversal` | `PinchoutTailError` | 6 | 0.026391 | 0.021083 (PASS) | **0.020747** | 0.029030 | **PASS** |
| `VelocityReversal` | `AVOInterceptError` | 12 | 0.464177 | 0.493383 (FAIL) | **0.421776** | 0.487386 | **PASS** |
| `VelocityReversal` | `AVOGradientError` | 12 | 0.683578 | 0.715841 (PASS) | **0.623769** | 0.717757 | **PASS** |

### 4.2 跨 9 模型 246 对全集中的非完美项、触边与目标失配诊断 (Transparent Disclosure of Edge Cases)

- **历史主基准 (`126` 对)**：`CombinedUncertainty` 在全部 **43/43** 项工程门槛上 **100% 通过 (`0` failed checks)**，尺度网格触边次数 (`ScaleEndpoint != 0`) 为 **`0/126`**。
- **扩展 9 模型多校准全集 (`246` 对，115 项检查)**：在每次校准种子均真实重建算子后，`CombinedUncertainty` 在全部 15 个 `Suite/Model` 组别的平均 `RelativeError` 上均取得正改善（`+3.08%` 至 `+16.10%`），通过 **114/115** 项检查（原先 `MultiCalMain/SoftLowQ/AVOGradientError` 已完全通过）。唯一未通过的 1 项局部检查如实披露如下：
  - `NewGeologyMain/Holdout_ThinInterbed/WeakEventError`：`{"baseline": 1.4179925427020568, "candidate": 1.5124526792886868, "allowance": 0.07089962713510284, "n": 3}`（盲测薄互层弱反射窗局部波动；同一模型下 `NoiseUncertaintyOnly` 100% 通过该局部检查）。

## 5. 四项显式模型误差与有色噪声压力测试集 (Explicit Model-Error & Colored-Noise Stress Suites)

其中 `StressVelocityError` 采用真实 $+4\%$ 速度模型扰动 (`vp_st = 1.04 * vp`, `h_st = 1.035 * h`) 并重新执行 `(Q, fpeak)` 校准与全套反演算子重建：

| 压力测试集 (Stress Suite) | 物理失配机制 (Mismatch Mechanism) | Frozen 均值 | StrongerDamping (η=1.15) | NoiseUncertaintyOnly (P1-A) | QModelErrorOnly (P1-B) | CombinedUncertainty (P1-C) |
|---|---|---:|---:|---:|---:|---:|
| `StressNonRickerPhase` | 真实子波含 $30^\circ$ 常相位旋转（校准仍假设零相位 Ricker） | 0.8495 | 0.8461 (+0.40%) | 0.7792 (+8.27%) | 0.8427 (+0.79%) | **0.7783 (+8.38%)** |
| `StressVelocityError` | 公开速度模型存在 $+4\%$ 偏差并重新校准 $(Q, f_{\text{peak}})$ 与重建传播算子 | 0.9971 | 0.9606 (+3.66%) | 0.8475 (+15.00%) | 1.0014 (-0.43%) | **0.8767 (+12.08%)** |
| `StressResidualNMO` | 远偏移距存在最高 $4\text{ ms}$ 抛物线残余动校正拉伸 | 0.8000 | 0.8023 (-0.29%) | 0.7293 (+8.84%) | 0.7976 (+0.29%) | **0.7278 (+9.02%)** |
| `StressColoredNoise` | 观测背景为 AR(1) 时间相关有色噪声 ($\rho=0.55$) | 0.7482 | 0.7573 (-1.22%) | 0.6974 (+6.79%) | 0.7246 (+3.15%) | **0.6916 (+7.56%)** |

## 6. 多次独立校准重复性与真值隔离/同先验公平性审计 (Calibration Repeatability & Fairness Audits)

| 模型组别 (Group) | 地质模型 (Model) | EffectiveQ (均值 ± 标准差) | EstimatedFrequency (Hz) | 校准对数Q不确定度 `q_uncertainty` | 校准对数频率不确定度 `f_uncertainty` | 权衡相关系数 `log_qf_corr` |
|---|---|---|---|---:|---:|---:|
| `historical` | `SoftLowQ` | 27.36 ± 0.31 | 39.37 ± 0.37 | 0.1114 ± 0.0070 | 0.1399 ± 0.0069 | -0.7380 |
| `historical` | `VelocityReversal` | 67.89 ± 0.16 | 76.76 ± 0.25 | 0.1180 ± 0.0023 | 0.0887 ± 0.0004 | -0.7626 |
| `historical` | `ThinHighFrequency` | 49.91 ± 0.68 | 65.17 ± 0.61 | 0.1333 ± 0.0031 | 0.1707 ± 0.0052 | -0.7119 |
| `development` | `Dev_ModerateQ` | 44.46 ± 0.28 | 47.55 ± 0.29 | 0.1390 ± 0.0011 | 0.1303 ± 0.0012 | -0.7676 |
| `development` | `Dev_DeepContrast` | 50.49 ± 0.45 | 66.33 ± 0.39 | 0.1091 ± 0.0052 | 0.1033 ± 0.0034 | -0.7627 |
| `heldout_new_geology` | `Holdout_StrongAtten` | 26.70 ± 0.10 | 44.54 ± 0.32 | 0.0927 ± 0.0062 | 0.1541 ± 0.0057 | -0.6379 |
| `heldout_new_geology` | `Holdout_HighContrastReversal` | 63.37 ± 0.19 | 79.47 ± 0.35 | 0.0918 ± 0.0016 | 0.0858 ± 0.0013 | -0.6301 |
| `heldout_new_geology` | `Holdout_ThinInterbed` | 47.00 ± 0.33 | 67.32 ± 0.49 | 0.1264 ± 0.0012 | 0.1507 ± 0.0044 | -0.7947 |
| `heldout_new_geology` | `Holdout_DeepLowQ` | 53.52 ± 0.49 | 63.52 ± 0.43 | 0.1951 ± 0.0036 | 0.1441 ± 0.0010 | -0.8254 |

- **真值隔离审计 (`results/truth_isolation_audit.json`)**：静态扫描 `core_p1/` 全部 `20` 个反演与校准 MATLAB 源码文件及 Python 双生求解器，公开输入严格限制在白名单，零真值访问；动态投毒测试反演差异 `0.0`。
- **Direct 与 Split 同先验/同预算公平性审计 (`results/matched_prior_fairness_audit.json`)**：跨全部 `1722` 对 `(Direct, Split)` 配对运行，两路线共享完全相同的 `EffectiveQ`、`QUncertainty`、`FUncertainty`、`sigma`、先验算子 `H`、模态权重 $R(s)$ 与尺度因子 `Scale`（最大差 `0.0`），且所有候选与无结构分支均满足 `FinalFlag == 0` 与 `NormalResidual <= 2e-7`。
- **Windows 本机共享 `ResearchMATLAB` 一键复现命令**：
  ```powershell
  python scripts/run_matlab_p1.py --action check --output results/matlab_live_check
  python scripts/run_matlab_p1.py --action validate --output results/matlab_live_validate
  python scripts/run_matlab_p1.py --action p1_experiment --output results/matlab_live_full
  ```
