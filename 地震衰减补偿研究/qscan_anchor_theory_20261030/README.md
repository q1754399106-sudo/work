# 子波锚定理论与锚点实验（2026-10-30）

理论推导 `THEORY_ZH.md`；开跑前写下的预测 `PLAN_ZH.md`（提交 1ac6e24）。

- 结论：QScan 选出的 1/c = (1 − κ)/c_true + κ/c_a + β。κ = a₀⟨a⟩/⟨a²⟩ 由模型 Q 和分析窗算得 0.675；锚点实验（实验 12 的 10 个剖面 × 锚点 0.7/1.0/1.4 × λ 1/2、1/16）测得 K = 0.66、0.71，A + K = 1.00、0.95，20/20 单调。现流程斜率 1 − κ ≈ 0.33（实测 0.19–0.30），权重只平移 β，所以调权重救不了压缩。
- 锚窗检验（预测先提交 91dee11）：把估子波的窗口移到 1.8–2.2 s、2.4–2.8 s，预测 κ = 0.96、1.14，实测 K = 0.98、1.12；最深窗口下选出的 c 与真 c 反向相关，与预测一致。四条判据全部满足；A 比 1 − κ 低 0.09–0.18。
- 补充（预测先提交 50470f2）：QScan-W 的目标函数曲率只有固定锚点的约 1/8（理论 8.0，实测中位数 9.8、6.5），与它偏差随数据漂移一致；其 1/c 斜率 1.2–1.4 的具体来源未单独证实。子波归一化的猜测在代码上被排除。
- 噪声检验（预测先提交 26d3ee7）：120 次扫描、噪声 ×0.25–×2。κ 与噪声无关（深窗 K 1.16–1.25，默认 0.65–0.69）；深锚窗 A + K 偏低主要来自噪声（低噪声 0.96–0.98，×2 时 0.82）；但 QScan-W 斜率高于 1 **不是噪声造成的**（×0.25、×0.5 选择相同，斜率都是 1.28），预测 1 不满足，其余三条满足。事后观察：低噪声下 QScan-W 近似常数比例偏差 c_sel ≈ 0.7 c_true（log 斜率 0.94）。
- 事后检验“偏差是否乘性”（`scripts/fit_ratio.py`，`analysis/fit_ratio.json`）：不成立。三组数据都是“共用斜率、λ 只平移”的加性模型拟合更好（自助 99–100%）；log 斜率随 λ 在 0.84–1.58 间变，只在截距过零时接近 1。稳定的是 1/c 斜率 A ≈ 1.2–1.4，与 λ 无关，所以超出来自失配项而非稀疏权重，具体来源未定。
- 截短分析窗检验（预测先提交 1433d09；60 次扫描，新种子 9941–9944）：深层来源被否定，分析窗从 1.2–3.0 s 截到 1.2–2.2 s 后 QScan-W 斜率 1.07→1.28，不降反升；对照斜率 1.07，未复现 1.1–1.5；截短窗下预测 κ = 0.804，实测 K = 0.81（κ 公式第五次通过）。事后观察：超出随曲率降低而变大，大小随数据批次在 0.07–0.39 间变，更像是数据相关的小扰动被低曲率放大。

- 迭代反馈检验（2026-11-01，预测先提交 a7fb41e）：把选出的 Q 作为下一步的锚点，24 个数据集 8 步全部单调收敛，增量比 0.689（理论 κ = 0.675）；8 步后与一步反解 (v_1 − κ)/(1 − κ) 的中位差 0.012；截距 0.107 > 0（理论 β/(1 − κ)），精度不比第一步好。P1、P4、P5 满足，P3（斜率 0.8–1.2）不满足：实测 0.76，是第一步斜率 0.24 被放大 1/(1 − κ) 倍的结果。mid / deep 组（κ = 0.964、1.14，各 6 个数据集）预登记的 P2 未通过：离开线性区（v ≳ 1.7）后迭代停滞或振荡；事后只看线性区，增量比 0.99、1.17，与 κ 相符（`analysis/iter_posthoc.json`）。
- 剖面形状检验（2026-11-01，预测先提交 1951185）：上覆层强吸收剖面 S 上 κ 预测 0.880、实测 0.856，现流程斜率 0.14（预测 0.12）；深层强吸收剖面 D 上预测 0.294、实测 0.598，**不满足**。截短分析窗的事后探索显示 D 上 2.2 s 以下贡献很小（截到 2.2 s 时实测 0.80、均匀权重 0.75）。κ 公式的适用边界：深层强吸收时低估压缩。
- 两层参数化的探索（开发种子）：深层 Q 的 lasso 目标随 v₂ 线性下降到 v₂ = 0，没有内部极小，深层 Q 总选"不衰减"。
- 论文初稿：`PAPER_DRAFT_ZH.md`。全部事先预测与实测：`figures/kappa_all.png`。

- 实测检验（2026-11-01，预测先提交 9cfb36f）：两个锚点各扫描一次即可在没有真值的数据上测 κ̂。Volve 10 条实测测线 κ̂ = 0.678（预测 0.675，测线间 0.669–0.690），深锚窗 1.19（预测 1.14），F0–F2 全部满足。κ̂ 已知后，实测上 QScan 与 SR 的差距主要来自 β（约 0.29）。
- 两锚点诊断检验（预测先提交 d6deecf）：剖面 D2（深层较强吸收）上公式 0.489、实测 0.624，剖面 M（中部吸收）上 0.595 对 0.614；两个剖面上 1 − κ̂ 都预测了压缩斜率。G1–G3 全部满足。公式只在衰减集中于深层时失效。

- 全工区实测（预测先提交 ca44394）：26 条 SR 测线上 κ̂ = 0.681（26/26 在 0.657–0.705），β̂ = 0.28（26/26 为正，测线间近似恒定），K1、K2、B1 全部满足；c_SR 在 0.85–1.37 之间变化时，现流程 QScan 只在 0.74–0.83 之间变化。

| 路径 | 内容 |
|---|---|
| `scripts/fit_existing.py` | 用实验 12、14 已有结果在 1/c 空间拟合（事后） |
| `scripts/anchor.py` | 锚点实验（固定锚点估子波的 QScan） |
| `scripts/curvature.py` | 曲率比检验（`PLAN_CURV_ZH.md`），`analysis/curvature.json` |
| `scripts/anchor2.py`、`scripts/analyze_window.py`、`scripts/fig_window.py` | 锚窗检验（`PLAN_WINDOW_ZH.md`），`results2/`、`analysis/window.json`、`figures/window_kappa.png` |
| `scripts/noise.py`、`scripts/analyze_noise.py` | 噪声检验（`PLAN_NOISE_ZH.md`），`results3/`、`analysis/noise.json`、`logs/` |
| `scripts/trunc.py`、`scripts/analyze_trunc.py` | 截短分析窗检验（`PLAN_TRUNC_ZH.md`），`results4/`、`analysis/trunc.json`、`logs/trunc_*` |
| `scripts/fit_ratio.py` | 事后检验 QScan-W 偏差是否乘性，`analysis/fit_ratio.json` |
| `scripts/analyze_anchor.py` | 锚点实验拟合、自助区间、单调性 |
| `scripts/fig_theory.py` | `figures/anchor_theory.png` |
| `scripts/env15.py` | 导入实验 7–13 的程序（同实验 14） |
| `results/` | 每次扫描一个 JSON（选择、目标值、收敛） |
| `analysis/` | `fit_existing.json`、`anchor.json` |
| `scripts/iterate.py`、`scripts/analyze_iter.py`、`scripts/fig_iter.py` | 迭代反馈检验（`PLAN_ITER_ZH.md`），`results5/`、`analysis/iter.json`、`figures/iter.png`、`logs/iter_*` |
| `scripts/qprofile.py`、`scripts/analyze_profile.py` | 剖面形状检验（`PLAN_PROFILE_ZH.md`），`results7/`、`analysis/profile.json`、`logs/prof_*` |
| `scripts/dtrunc.py` | D 剖面截短分析窗（事后探索），`results8/`、`analysis/dtrunc.json` |
| `scripts/layered.py` | 两层参数化扫描（冒烟后发现深层不可辨，未正式运行；改用 `probe_layer.py` 沿深层单轴扫描） |
| `scripts/fig_kappa_all.py` | `figures/kappa_all.png`：全部 κ 预测与实测 |
| `PAPER_DRAFT_ZH.md` | 论文初稿 |
| `scripts/probe_layer.py`、`scripts/probe_active.py` | 两层参数化与非零系数加权的事后探索，`logs/probe_layer.log`、`logs/probe_active.log` |
| `scripts/fig_schematic.py` | `figures/schematic.png`：论文图 1（吸收剖面、窗口与 κ 随子波窗深度的变化） |
| `scripts/posthoc_iter.py` | 迭代检验 mid / deep 组只用可靠步的增量比（事后），`analysis/iter_posthoc.json` |
| `scripts/field_kappa.py`、`scripts/analyze_field.py` | 实测测线上的两锚点 κ̂（`PLAN_FIELD_ZH.md`），`results9/`、`analysis/field.json`、`logs/field_*` |
| `scripts/qprofile2.py`、`scripts/analyze_profile2.py` | 两锚点诊断与剖面 D2、M（`PLAN_PROFILE2_ZH.md`），`results10/`、`analysis/profile2.json`、`logs/prof2_*` |
| `scripts/field_kappa2.py`、`scripts/analyze_field2.py`、`scripts/fig_field.py` | 全工区 26 条 SR 测线的 κ̂ 与 β（`PLAN_FIELD2_ZH.md`），`results11/`、`analysis/field2.json`、`figures/field.png`、`logs/field2_*` |
