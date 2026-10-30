# 子波锚定理论与锚点实验（2026-10-30）

理论推导 `THEORY_ZH.md`；开跑前写下的预测 `PLAN_ZH.md`（提交 1ac6e24）。

- 结论：QScan 选出的 1/c = (1 − κ)/c_true + κ/c_a + β。κ = a₀⟨a⟩/⟨a²⟩ 由模型 Q 和分析窗算得 0.675；锚点实验（实验 12 的 10 个剖面 × 锚点 0.7/1.0/1.4 × λ 1/2、1/16）测得 K = 0.66、0.71，A + K = 1.00、0.95，20/20 单调。现流程斜率 1 − κ ≈ 0.33（实测 0.19–0.30），权重只平移 β，所以调权重救不了压缩。
- 未解释：每个 c 重估子波（QScan-W）时 1/c 空间斜率 1.2–1.4，而不是 1。

| 路径 | 内容 |
|---|---|
| `scripts/fit_existing.py` | 用实验 12、14 已有结果在 1/c 空间拟合（事后） |
| `scripts/anchor.py` | 锚点实验（固定锚点估子波的 QScan） |
| `scripts/analyze_anchor.py` | 锚点实验拟合、自助区间、单调性 |
| `scripts/fig_theory.py` | `figures/anchor_theory.png` |
| `scripts/env15.py` | 导入实验 7–13 的程序（同实验 14） |
| `results/` | 每次扫描一个 JSON（选择、目标值、收敛） |
| `analysis/` | `fit_existing.json`、`anchor.json` |
