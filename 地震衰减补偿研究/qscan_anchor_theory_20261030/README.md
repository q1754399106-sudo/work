# 子波锚定理论与锚点实验（2026-10-30）

理论推导 `THEORY_ZH.md`；开跑前写下的预测 `PLAN_ZH.md`（提交 1ac6e24）。

- 结论：QScan 选出的 1/c = (1 − κ)/c_true + κ/c_a + β。κ = a₀⟨a⟩/⟨a²⟩ 由模型 Q 和分析窗算得 0.675；锚点实验（实验 12 的 10 个剖面 × 锚点 0.7/1.0/1.4 × λ 1/2、1/16）测得 K = 0.66、0.71，A + K = 1.00、0.95，20/20 单调。现流程斜率 1 − κ ≈ 0.33（实测 0.19–0.30），权重只平移 β，所以调权重救不了压缩。
- 锚窗检验（预测先提交 91dee11）：把估子波的窗口移到 1.8–2.2 s、2.4–2.8 s，预测 κ = 0.96、1.14，实测 K = 0.98、1.12；最深窗口下选出的 c 与真 c 反向相关，与预测一致。四条判据全部满足；A 比 1 − κ 低 0.09–0.18。
- 补充（预测先提交 50470f2）：QScan-W 的目标函数曲率只有固定锚点的约 1/8（理论 8.0，实测中位数 9.8、6.5），与它偏差随数据漂移一致；其 1/c 斜率 1.2–1.4 的具体来源未单独证实。子波归一化的猜测在代码上被排除。
- 噪声检验（预测先提交 26d3ee7）：120 次扫描、噪声 ×0.25–×2。κ 与噪声无关（深窗 K 1.16–1.25，默认 0.65–0.69）；深锚窗 A + K 偏低主要来自噪声（低噪声 0.96–0.98，×2 时 0.82）；但 QScan-W 斜率高于 1 **不是噪声造成的**（×0.25、×0.5 选择相同，斜率都是 1.28），预测 1 不满足，其余三条满足。事后观察：低噪声下 QScan-W 近似常数比例偏差 c_sel ≈ 0.7 c_true（log 斜率 0.94）。
- 事后检验“偏差是否乘性”（`scripts/fit_ratio.py`，`analysis/fit_ratio.json`）：不成立。三组数据都是“共用斜率、λ 只平移”的加性模型拟合更好（自助 99–100%）；log 斜率随 λ 在 0.84–1.58 间变，只在截距过零时接近 1。稳定的是 1/c 斜率 A ≈ 1.2–1.4，与 λ 无关，所以超出来自失配项而非稀疏权重，具体来源未定。
- 截短分析窗检验（预测先提交 1433d09；60 次扫描，新种子 9941–9944）：深层来源被否定，分析窗从 1.2–3.0 s 截到 1.2–2.2 s 后 QScan-W 斜率 1.07→1.28，不降反升；对照斜率 1.07，未复现 1.1–1.5；截短窗下预测 κ = 0.804，实测 K = 0.81（κ 公式第五次通过）。事后观察：超出随曲率降低而变大，大小随数据批次在 0.07–0.39 间变，更像是数据相关的小扰动被低曲率放大。

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
