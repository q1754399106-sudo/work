"""Relative combination error of every criterion against QScan + Standard (A0), per scene, and the
selection error in grid steps (analysis/criteria_scores.csv)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
s = pd.read_csv(ROOT / "analysis" / "criteria_scores.csv")
crit = ["A1", "A2", "A3", "C1", "C2", "C3", "C4"]
Q0, Q07, LOW, INK, MUTED, SURFACE, GRID = "#2a78d6", "#eb6834", "#a3a29d", "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"

plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={"width_ratios": [3, 2]}, facecolor=SURFACE)
rng = np.random.default_rng(0)
for i, c in enumerate(crit):
    for tag, color, label in (("Q100_14dB", Q0, "Q0, 14 dB"), ("Q070_14dB", Q07, "Q0 x 0.7, 14 dB"), ("Q100_08dB", LOW, "Q0, 8 dB (descriptive)")):
        v = 100 * s[(s.criterion == c) & (s.tag == tag)].rel_vs_A0.to_numpy()
        ax.plot(i + rng.uniform(-0.18, 0.18, v.size), np.clip(v, -80, 150), "o", ms=3.5, color=color, alpha=0.85,
                label=label if i == 0 else None)
    p = s[(s.criterion == c) & s.primary]
    m = 100 * p.rel_vs_A0.mean()
    ax.plot([i - 0.28, i + 0.28], [min(m, 150)] * 2, color=INK, lw=1.6, label="mean of 20 primary scenes" if i == 0 else None)
    ax.text(i, -76, f"mean {m:+.0f}%", ha="center", va="bottom", color=INK, fontsize=7)
ax.axhline(0, color=MUTED, lw=0.8)
ax.axhspan(-80, -10, color=Q0, alpha=0.05, lw=0)
ax.set_xticks(range(len(crit)), crit)
ax.set_ylim(-80, 150)
ax.set_ylabel("combination MSE change vs QScan (A0), %")
ax.set_title("Each criterion's choice, scored after coherence Standard (values beyond -80/+150 clipped)", color=INK, fontsize=9, loc="left")
ax.legend(frameon=False, fontsize=7, loc="upper left", ncol=2)

for i, c in enumerate(["A0"] + crit):
    p = s[(s.criterion == c) & s.primary]
    bx.plot(i + rng.uniform(-0.18, 0.18, len(p)), p.mean_step_error, "o", ms=3.5, color=Q0 if c != "A0" else INK, alpha=0.8)
bx.axhline(0, color=MUTED, lw=0.8)
bx.set_xticks(range(len(crit) + 1), ["A0"] + crit)
bx.set_ylabel("mean selection error (grid steps, + = Q too high)")
bx.set_title("Selection error on the 20 primary scenes", color=INK, fontsize=9, loc="left")
for a in (ax, bx):
    a.set_facecolor(SURFACE)
    a.grid(axis="y", color=GRID, linewidth=0.5)
    a.set_axisbelow(True)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
fig.text(0.01, 0.01, "F1, F2, G1-G8; replicates 1-6. Shaded band: the W1 mean-improvement target (<= -10%).", color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"criteria_combo.{ext}", dpi=200, facecolor=SURFACE)
print("written", FIG / "criteria_combo.png")
