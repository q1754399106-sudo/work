# `uq_rtm_closure_20261125` 预登记冻结执行方案（`WP0–WP5`，门控 `R0–R7`）

- **日期**：2026-11-26
- **分支与隔离工作树**：`agy/uq-rtm-closure-20261125`（`D:\ChatGPT\研究\.claude\worktrees\uq-rtm-closure-20261125`，自 `origin/main` `09710f3` 切出）
- **任务书依据**：`reports/TASK_SPEC_UQ_RTM_REMEDIATION_AND_CLOSURE_20261125_ZH.md`
- **只读上游数据**：`uq_rtm_20261117/cache/fld`、`study-research-docs-plan-c4d4c6/.../uq_rtm_20261117/cache/fld45`、`D:\tmp_uqrtm_unified_20261124\`、`uq_rtm_unified_20261124/`（历史证据原封不动）

---

## 1. 冻结前已完成前置工作包

1. **`WP0`（科研诚信机械审计门神）**：
   - 脚本：`scripts/audit_report_integrity.py`、`scripts/render_full_tables.py`
   - 门神自测（`logs/wp0_audit_selftest.log`）：针对旧报告 `uq_rtm_unified_20261124/RESULTS_ZH.md` 自动捕获 **38 项违规**（含 `F-14` 漏列、`U4b` 漏法、`U1b/U3b` 漏 `SA` 行、4–45 Hz `G_HF=1.054<1.10` 休眠未披露与微差夸大词），退出码 `1`。
2. **`WP1`（零重算 / 轻量异常根因诊断 `D1–D4`）**：
   - 脚本：`scripts/diag_anomalies_wp1.m`；产出：`results/wp1_diagnostics.json`、`figures/fig_wp1_diagnostics.png`、`logs/wp1_diagnostics.log`。
   - 四项根因全部定量闭合（`D1` `F-14` `-55~-60 m` 深度滞后被 `±30 m` 窗截断在负旁瓣；`D2` `RAW>DSA` 中 `21–23%` 来自 `PSAE` 斜率差、`77–79%` 来自 `F-5` 覆次减半与上覆画弧 `63.4% → 3.9%` 权衡；`D3` `25 Hz` 锁定 `-20 m` 与 `RAW_R0` 一致而 `100 Hz` 借人工相位漂移偶然过零；`D4` `U4b` 中 `DSA_UQW1b_PAnc_25` 对 2D 弹性真值 `E1_vs_el` 最高 `0.8424`，`well_cc` 倒挂源于无振幅权重 1D 反射率相关 + `B=300` 在 20–25 Hz 略偏保守）。
3. **`WP4` 预冻结 SNR 冒烟校验（`logs/wp4_snr_smoke.log`）**：
   - 经 10 频代表性频网实测校准，锁定 `U6/U8/U9: snr_db = +0.5 dB`、`U7: snr_db = -3.0 dB`，对应全 34 频 `RAW_R0 W_deep SNR_img ≈ +4.1 ~ +4.7 dB`（全部稳落 `[+2.0, +6.0] dB` 区间中心，杜绝 `-14 dB` 定标失误）。

---

## 2. 预登记待执行生产工作包（`WP2–WP4`）

1. **`WP3`（P0 `PSAE` 归因对照，Volve 4–25 Hz，85 频）**：
   - 脚本：`scripts/run_closure_volve25_psae_ctrl.m` → `results/volve25_psae_ctrl.json`
   - 15 组全量对照（路线 A 未锚定、路线 B1/B2 去 `(2πf)^2` 零 `PSAE`、路线 C `PSAE`、路线 D 6 Hz 平滑浅层谱反褶积，及 `PH`/`IQ70`/`W1b` 公平竞争者）。
2. **`WP4`（P1 全新独立合成套件 `U6–U9`，种子 `20261401–20261404`，34 频 `4:1:37 Hz`）**：
   - 脚本：`scripts/run_closure_synth_u6_u9.m` → `results/synth_U6.json`–`U9.json`
   - 涵盖气烟囱缓倾（`U6`）、22° 陡倾逆断层（`U7`）、空间错位+尺度失配非标量 Q 误差（`U8`，`x=970 m, rx×1.4, m∈{0.6,1.0,1.4}`）、深浅反向非标量 Q 误差（`U9`，浅 `1.5Q`/深 `0.6Q`）。
3. **`WP2`（P0 Volve 4–45 Hz 真无混叠同域 `INV3_DSA` 合流重跑，全 165 频 `df=0.25 Hz`）**：
   - 脚本：`scripts/run_closure_volve45.m` → `results/volve45_closure.json`
   - 两遍轻量架构（`D:\tmp_uqrtm_closure_20261125\` 峰值 `< 2.5 GB`）：Pass 1 提取并落盘 `INV3_DSA`、`INV3_RAW_FULL`、`INV3_DSA_NOW2`；Pass 2 真实接线注入全部 18 组对照（含旧混叠 `S4Hiv` 复现对照、真同域 `TrueDSAiv`、去 `(2πf)^2` `NoW2`、浅层 6 Hz 反褶积 `Decon`、`SameWeq` 等）。

---

## 3. 预登记冻结机械门控（`R0–R7`，`scripts/judge_closure.py` 冻结不改）

- **`R0_Integrity_Gate`**：`audit_report_integrity.py` 在旧报告上报错（已过）且在新报告 `RESULTS_ZH.md` 上退出码 `0`（一票否决）。
- **`R1_Active_Broadband_Operator`**：WP2 中 `DSA_UQW1b_TrueDSAiv_PAnc_25 (m=1)` 在 `28–40 Hz` 相对 `DSA_R0_Eq` 中位振幅增益比 $G_{\text{HF}} \ge 1.25$。
- **`R2_Volve45_TrueConfluence_Closure`**：① `Smile_Frac <= 0.10` 且相对 Stage 4H 压降 $\ge 85\%$；② `F5_fix` 5-$m$ 均分 $\ge \text{RAW\_R0} + 0.040$ 且 5/5 $m > \text{RAW\_R0}$；③ 储层 `[6–38 Hz]` 谱斜率 $|\beta_{\text{res}}|$ 较 `DSA_R0_Eq` 至少改善 30%。
- **`R3_Volve45_TrueUQ_vs_Det_At_m06`**：$m=0.6$ 下 `DSA_UQW1b_TrueDSAiv_PAnc_25` 严格高于 `DSA_IQ70_PAnc_25`（$\ge +0.020$）且 $\ge$ 同域真 `INV3_DSA` 的 `DSA_W1b_TrueDSAiv_PAnc_25`。
- **`R4_PSAE_vs_NoW2_Attribution`**：WP3（4–25 Hz）与 WP2（4–45 Hz）完整量化“去 `(2πf)^2` 零 `PSAE`”与“`(2πf)^2 + PSAE`”在 `F5`、`F14`、`Res_Slope` 上的精确差异与归因占比。
- **`R5_Synth_U6U7_DualMetric_And_SA_Tradeoff`**：`U6–U7` 上① 定理 3 `DSA B=300` 2/2 落入 `±1.0` 档；② `E_amp` 较 `DSA_PH_25` 降低 $\ge 30\%$（2/2）；③ $m=0.5$ 的 `E1` 较 `max(DSA_IQ70, DSA_W1b)` 高 $\ge +0.030$（2/2）；④ 全量如实报告 `SA` vs `DSA`。
- **`R6_NonScalar_Spatial_Q_Error_U8U9`**：`U8–U9` 上 `DSA_UQW1b_25` 的 `E1` 均分严格高于 `DSA_IQ70_25` 且最差工况 $\ge \text{RAW\_R0} - 0.020$。
- **`R7_Anomaly_Resolution`**：`WP1` 的 `D1–D4` 定量归因 JSON 与诊断图完整落盘。
