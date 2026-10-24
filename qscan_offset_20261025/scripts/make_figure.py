"""Error curves of the eight development models relative to their own optimum, and the
relative change of the shifted selections against QScan + Standard (analysis/offset_*.csv)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
c = pd.read_csv(ROOT / "analysis" / "offset_curves.csv")
s = pd.read_csv(ROOT / "analysis" / "offset_summary.csv")
scenes = [("Q100_14dB", "Q0, 14 dB"), ("Q100_08dB", "Q0, 8 dB"), ("Q070_14dB", "Q0 x 0.7, 14 dB")]
LINE, QSCAN, OPT, INK, MUTED, SURFACE = "#2a78d6", "#c0392b", "#0b0b0b", "#0b0b0b", "#52514e", "#fcfcfb"
SHIFT = {"rel_S1": ("#7aa9e0", "shift 1 step"), "rel_S2": ("#2a78d6", "shift 2 steps"), "rel_opt": ("#a3a29d", "best fixed factor")}

plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
fig, axes = plt.subplots(1, 4, figsize=(12, 3.6), gridspec_kw={"width_ratios": [3, 3, 3, 3.4]}, facecolor=SURFACE)
for ax, (tag, title) in zip(axes, scenes):
    ax.set_facecolor(SURFACE)
    for model in sorted(c.model.unique()):
        k = c[(c.model == model) & (c.scene == tag)].sort_values("q_ratio_to_Q0")
        row = s[(s.model == model) & (s.scene == tag)].iloc[0]
        rel = k.mse / k.mse.min()
        ax.plot(k.q_ratio_to_Q0, rel, color=LINE, lw=1, alpha=0.55)
        at = k[np.isclose(k.q_ratio_to_Q0, row.q_qscan_mode)]
        ax.plot(at.q_ratio_to_Q0, at.mse / k.mse.min(), "o", color=QSCAN, ms=4)
    ax.axvline(1, color=MUTED, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xticks([0.5, 0.7, 1, 1.4, 2], ["0.5", "0.7", "1", "1.4", "2"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(0.55 if tag != "Q070_14dB" else 0.4, 2.3 if tag != "Q070_14dB" else 1.6)
    ax.set_ylim(0.95, 2.0)
    ax.set_xlabel("effective Q / upscaled Q0")
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    ax.grid(color="#e4e3df", linewidth=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
axes[0].set_ylabel("MSE / model's minimum MSE")
axes[0].plot([], [], color=LINE, label="one model (G1-G8)")
axes[0].plot([], [], "o", color=QSCAN, ms=4, label="QScan's modal choice")
axes[0].legend(frameon=False, loc="upper center")

ax = axes[3]
ax.set_facecolor(SURFACE)
x = np.arange(len(scenes))
for j, (col, (color, label)) in enumerate(SHIFT.items()):
    for i, (tag, _) in enumerate(scenes):
        v = 100 * s[s.scene == tag][col].to_numpy()
        ax.plot(np.full(v.size, i + (j - 1) * 0.25), v, "o", color=color, ms=3.5, alpha=0.85, label=label if i == 0 else None)
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xticks(x, ["Q0 14 dB", "Q0 8 dB", "Q0x0.7 14 dB"])
ax.set_ylabel("MSE change vs QScan + Standard (%)")
ax.set_title("Moving QScan's choice down", color=INK, fontsize=9, loc="left")
ax.grid(axis="y", color="#e4e3df", linewidth=0.5)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(frameon=False, loc="lower left", fontsize=7)
fig.text(0.01, 0.01, "Each curve: Standard MSE (R = 6) over the 14 factors, divided by that model's minimum; curves above 2.0 are clipped.",
         color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"offset_curves.{ext}", dpi=200, facecolor=SURFACE)
print("written", FIG / "offset_curves.png")
