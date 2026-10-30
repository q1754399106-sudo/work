# 子波锚定理论与锚点实验（2026-10-30）

理论推导 `THEORY_ZH.md`；开跑前写下的预测 `PLAN_ZH.md`（提交 1ac6e24）。

- 结论：QScan 选出的 1/c = (1 − κ)/c_true + κ/c_a + β。κ = a₀⟨a⟩/⟨a²⟩ 由模型 Q 和分析窗算得 0.675；锚点实验（实验 12 的 10 个剖面 × 锚点 0.7/1.0/1.4 × λ 1/2、1/16）测得 K = 0.66、0.71，A + K = 1.00、0.95，20/20 单调。现流程斜率 1 − κ ≈ 0.33（实测 0.19–0.30），权重只平移 β，所以调权重救不了压缩。
- 锚窗检验（预测先提交 91dee11）：把估子波的窗口移到 1.8–2.2 s、2.4–2.8 s，预测 κ = 0.96、1.14，实测 K = 0.98、1.12；最深窗口下选出的 c 与真 c 反向相关，与预测一致。四条判据全部满足；A 比 1 − κ 低 0.09–0.18。
- 补充（预测先提交 50470f2）：QScan-W 的目标函数曲率只有固定锚点的约 1/8（理论 8.0，实测中位数 9.8、6.5），与它偏差随数据漂移一致；其 1/c 斜率 1.2–1.4 的具体来源未单独证实。子波归一化的猜测在代码上被排除。

| 路径 | 内容 |
|---|---|
| `scripts/fit_existing.py` | 用实验 12、14 已有结果在 1/c 空间拟合（事后） |
| `scripts/anchor.py` | 锚点实验（固定锚点估子波的 QScan） |
| `scripts/curvature.py` | 曲率比检验（`PLAN_CURV_ZH.md`），`analysis/curvature.json` |
| `scripts/anchor2.py`、`scripts/analyze_window.py`、`scripts/fig_window.py` | 锚窗检验（`PLAN_WINDOW_ZH.md`），`results2/`、`analysis/window.json`、`figures/window_kappa.png` |
| `scripts/analyze_anchor.py` | 锚点实验拟合、自助区间、单调性 |
| `scripts/fig_theory.py` | `figures/anchor_theory.png` |
| `scripts/env15.py` | 导入实验 7–13 的程序（同实验 14） |
| `results/` | 每次扫描一个 JSON（选择、目标值、收敛） |
| `analysis/` | `fit_existing.json`、`anchor.json` |
