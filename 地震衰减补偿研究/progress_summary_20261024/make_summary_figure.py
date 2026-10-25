"""Draw the progress-summary comparison figure from archived metric tables.

Reads the CSVs directly from git HEAD (works under sparse checkout) and writes
figures/method_comparison.png / .pdf next to this script. No experiment is rerun.

Run (repo root):
  python -I progress_summary_20261024/make_summary_figure.py
Needs NumPy, pandas, Matplotlib (attenuation_run_20261019/dependencies).
"""
import io
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "attenuation_run_20261019" / "dependencies"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

KEY = ["model", "condition", "Q_percent"]


def read_csv(path):
    local = ROOT / path
    if local.is_file():
        return pd.read_csv(local, encoding="utf-8-sig")
    raw = subprocess.run(["git", "show", f"HEAD:{path}"], cwd=ROOT,
                         capture_output=True, check=True).stdout
    return pd.read_csv(io.BytesIO(raw), encoding="utf-8-sig")


# 20 development scenes: four-layer (primary) and complex models, 2 noise levels, Q -30..+30%.
sr = read_csv("soft_refit_20261023/analysis/metrics.csv")
cd = read_csv("coherence_soft_20261023/analysis/development_metrics.csv")
bq = read_csv("bae_q_20261024/analysis/dev_metrics_all.csv")
dev = sr.pivot_table(index=KEY, columns="method", values="relative_error")
dev["coh_v1"] = cd[cd.method == "coh_1_8"].set_index(KEY).overall_error
bq30 = bq[bq.Q_percent.abs() <= 30].pivot_table(index=KEY, columns="method",
                                                values="relative_error")
dev["qscan"] = bq30.QScanRefit

# Unseen / sealed data from round 3 (F3/F4 generated after rules were frozen; V5 sealed complex model).
f34 = read_csv("round3_20261024/analysis/guard_metrics_F3_F4.csv")
v5 = read_csv("round3_20261024/analysis/guard_metrics_V5.csv")
# Original FixedSplit rerun on the same round3 observations (fixed_baseline_round3_20261024).
fx = read_csv("fixed_baseline_round3_20261024/analysis/fixed_metrics_F3_F4_V5.csv")
f34 = pd.concat([f34, fx[fx.set.isin(["F3", "F4"])]], ignore_index=True)
v5 = pd.concat([v5, fx[fx.set == "V5"]], ignore_index=True)

DEV_METHODS = [
    ("fixed", "最初分步法（10-19 基线）", "#7f7f7f"),
    ("unpruned", "FullCalRefit 未剪枝", "#1f77b4"),
    ("hard", "硬剪枝", "#9467bd"),
    ("soft", "软收缩", "#17becf"),
    ("coh_v1", "相干软收缩 v1", "#2ca02c"),
    ("qscan", "QScanRefit", "#ff7f0e"),
]
NEW_METHODS = [
    ("FixedSplit", "最初分步法（10-19 基线）", "#7f7f7f", "-"),
    ("FullCalRefit", "FullCalRefit", "#1f77b4", "-"),
    ("SoftGroupRefit", "软收缩", "#17becf", "-"),
    ("Standard", "相干软收缩标准版", "#2ca02c", "-"),
    ("QScanRefit", "QScanRefit", "#ff7f0e", "-"),
    ("QScan_Standard", "QScan + 相干标准版", "#d62728", "-"),
]

fig, axes = plt.subplots(2, 2, figsize=(16, 11.5))

# (a) Q correct, development scenes.
ax = axes[0, 0]
q0 = dev.xs(0, level="Q_percent")
cases = [("primary", "Moderate", "四层\n中等噪声"), ("primary", "Strong", "四层\n强噪声"),
         ("complex", "Moderate", "复杂\n中等噪声"), ("complex", "Strong", "复杂\n强噪声")]
x = np.arange(len(cases))
w = 0.8 / len(DEV_METHODS)
for i, (m, label, color) in enumerate(DEV_METHODS):
    vals = [q0.loc[(mo, co), m] for mo, co, _ in cases]
    bars = ax.bar(x + (i - (len(DEV_METHODS) - 1) / 2) * w, vals, w, label=label, color=color)
    if m in ("fixed", "coh_v1"):
        ax.bar_label(bars, fmt="%.2f", fontsize=7.5, padding=1)
ax.set_xticks(x, [c[2] for c in cases])
ax.set_ylabel("整体相对波形误差（越低越好）")
ax.set_title("(a) Q 正确时的精度：同一批 20 个开发场景", loc="left", fontsize=12)
ax.legend(fontsize=8.5, ncol=2, frameon=False)
ax.set_ylim(0, 1.18)
ax.grid(axis="y", alpha=0.3)

# (b) Q-error robustness, development scenes, mean of the four cases.
ax = axes[0, 1]
curve = dev.groupby(level="Q_percent").mean()
for m, label, color in DEV_METHODS:
    if m in ("hard", "soft"):
        continue
    ax.plot(curve.index, curve[m], "o-", color=color, label=label, lw=2)
ax.set_xlabel("所用 Q 相对真值 / 标定值的误差（%）")
ax.set_ylabel("4 个场景平均误差")
ax.set_title("(b) 对 Q 误差的稳健性：同一批开发场景", loc="left", fontsize=12)
ax.legend(fontsize=9, frameon=False)
ax.grid(alpha=0.3)
ax.set_xticks([-30, -10, 0, 10, 30])


def plot_new(ax, df, title):
    c = df.groupby(["Q_percent", "method"]).overall.mean().unstack()
    for m, label, color, ls in NEW_METHODS:
        ax.plot(c.index, c[m], marker="o", ls=ls, color=color, label=label, lw=2.2 if m == "QScan_Standard" else 1.6)
    ax.set_yscale("log")
    ax.set_xticks([-50, -30, -10, 0, 10, 30, 50])
    ax.set_xlabel("所用 Q 的误差（%）")
    ax.set_ylabel("平均整体误差（对数轴）")
    ax.set_title(title, loc="left", fontsize=12)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=9, frameon=False)


plot_new(axes[1, 0], f34, "(c) 冻结后新生成的独立正演 F3/F4（28 场景）")
plot_new(axes[1, 1], v5, "(d) 封存复杂模型 V5（14 场景；Q 0% 即地面资料标定 Q）")
for ax in axes[1]:
    ax.set_yticks([0.4, 0.5, 0.7, 1, 2, 4])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

fig.suptitle("地震 Q 补偿：方法演进与效果对比（合成数据，误差 = ‖输出 − 无吸收目标‖ / ‖目标‖）", fontsize=14)
fig.text(0.01, 0.005,
         "数据：soft_refit_20261023、coherence_soft_20261023（相干 v1 为开发集选参结果）、bae_q_20261024、round3_20261024。"
         "(c)(d) 中的最初分步法由 fixed_baseline_round3_20261024 在相同观测上补跑。QScanRefit 在 (a)(b) 中来自 bae_q 的同一批场景。",
         fontsize=8.5, color="#444")
fig.tight_layout(rect=(0, 0.02, 1, 0.96))
out = HERE / "figures"
out.mkdir(exist_ok=True)
fig.savefig(out / "method_comparison.png", dpi=200)
fig.savefig(out / "method_comparison.pdf")
print("wrote", out / "method_comparison.png")
