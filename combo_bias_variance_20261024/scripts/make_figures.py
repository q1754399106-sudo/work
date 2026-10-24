"""Bias^2 + variance of the four estimators in every scene (analysis/decomposition.csv)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
d = pd.read_csv(ROOT / "analysis" / "decomposition.csv")
d = d[d.region == "all"]
methods = [("FullCalRefit", "A"), ("Standard", "B"), ("QScanRefit", "C"), ("QScan_Standard", "D")]
panels = [("Q0, 14 / 8 dB", ["S01_F1_Q100_14dB", "S02_F2_Q100_14dB", "S03_F1_Q100_08dB", "S04_F2_Q100_08dB"]),
          ("Q0 x 0.7 (Q underestimated)", ["S05_F1_Q070_14dB", "S06_F2_Q070_14dB", "S07_F1_Q070_08dB", "S08_F2_Q070_08dB"]),
          ("Q0, 20 dB", ["S09_F1_Q100_20dB", "S10_F2_Q100_20dB"])]
BIAS, VAR, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), gridspec_kw={"width_ratios": [4, 4, 2]}, facecolor=SURFACE)
for ax, (title, scenes) in zip(axes, panels):
    ax.set_facecolor(SURFACE)
    ticks, labels = [], []
    for g, s in enumerate(scenes):
        for k, (m, short) in enumerate(methods):
            row = d[(d.scene == s) & (d.method == m)].iloc[0]
            x = g * 5 + k
            ax.bar(x, row.bias2, width=0.8, color=BIAS, edgecolor=SURFACE, linewidth=1)
            ax.bar(x, row["var"], bottom=row.bias2, width=0.8, color=VAR, edgecolor=SURFACE, linewidth=1)
            ticks.append(x)
            labels.append(short)
        model, snr = s.split("_")[1], s.split("_")[3]
        ax.text(g * 5 + 1.5, -0.12, f"{model} {snr}", transform=ax.get_xaxis_transform(), ha="center", va="top", color=INK)
    ax.set_xticks(ticks, labels)
    ax.tick_params(axis="x", length=0)
    ax.set_title(title, color=INK, fontsize=10, loc="left")
    ax.grid(axis="y", color="#e4e3df", linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
axes[0].set_ylabel("MSE / ||x*||^2")
handles = [plt.Rectangle((0, 0), 1, 1, color=BIAS), plt.Rectangle((0, 0), 1, 1, color=VAR)]
fig.legend(handles, ["bias^2", "variance"], loc="upper right", frameon=False, ncol=2)
fig.text(0.01, 0.01, "A FullCalRefit   B coherence Standard   C QScanRefit   D QScan + Standard   (R = 100 / 50 / 40 noise replicates)",
         color=MUTED, fontsize=8)
fig.tight_layout(rect=(0, 0.05, 1, 0.93))
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"bias_variance_2x2.{ext}", dpi=200, facecolor=SURFACE)
print("written", FIG / "bias_variance_2x2.png")
