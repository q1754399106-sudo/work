"""Figures of this round (descriptive; no judgment reads them).

figures/fd_curves.png/pdf  (a-c) Standard error curves of H1-H6 relative to each model's minimum,
                           with QScan's modal choice; (d) QScan offset and (e) noise-free optimum
                           at Q0/14 dB for the Kennett models G1-G8 and the 2-D FD models H1-H6.
figures/fd_gathers.png     one 2-D model: elastic target, noise-free attenuated data and their
                           difference (fault diffractions, gas-sand wedge).
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
c = pd.read_csv(ROOT / "analysis" / "fd_curves.csv")
s = pd.read_csv(ROOT / "analysis" / "fd_summary.csv")
g = pd.read_csv(ROOT / "cache" / "inputs" / "qscan_offset_20261025" / "analysis" / "offset_summary.csv")
scenes = [("Q100_14dB", "Q0, 14 dB"), ("Q100_08dB", "Q0, 8 dB"), ("Q070_14dB", "Q0 x 0.7, 14 dB")]
LINE, QSCAN, INK, MUTED, GRID, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})


def tidy(ax, title):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    ax.grid(color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


fig, axes = plt.subplots(1, 5, figsize=(14, 3.6), gridspec_kw={"width_ratios": [3, 3, 3, 2, 2]}, facecolor=SURFACE)
for ax, (tag, title) in zip(axes, scenes):
    for model in sorted(c.model.unique()):
        k = c[(c.model == model) & (c.scene == tag)].sort_values("q_ratio_to_Q0")
        row = s[(s.model == model) & (s.scene == tag)].iloc[0]
        ax.plot(k.q_ratio_to_Q0, k.mse / k.mse.min(), color=LINE, lw=1.2, alpha=0.6)
        at = k[np.isclose(k.q_ratio_to_Q0, row.q_qscan_mode)]
        ax.plot(at.q_ratio_to_Q0, at.mse / k.mse.min(), "o", color=QSCAN, ms=5, mec=SURFACE, mew=1)
    ax.axvline(1, color=MUTED, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xticks([0.5, 0.7, 1, 1.4, 2], ["0.5", "0.7", "1", "1.4", "2"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(0.55 if tag != "Q070_14dB" else 0.4, 2.3 if tag != "Q070_14dB" else 1.6)
    ax.set_ylim(0.95, 2.0)
    ax.set_xlabel("effective Q / upscaled Q0")
    tidy(ax, title)
axes[0].set_ylabel("MSE / model's minimum MSE")
axes[0].plot([], [], color=LINE, label="one 2-D model (H1-H6)")
axes[0].plot([], [], "o", color=QSCAN, ms=5, label="QScan's modal choice")
axes[0].legend(frameon=False, loc="upper center")

p = s[s.scene == "Q100_14dB"]
gp = g[g.scene == "Q100_14dB"]
for ax, col, title, ylab in [(axes[3], "offset_steps", "QScan offset, Q0 14 dB", "grid steps (+ = Q too high)"),
                             (axes[4], "q_best_noise_free", "Noise-free optimum, Q0 14 dB", "effective Q / upscaled Q0")]:
    for i, (frame, lab) in enumerate([(gp, "Kennett 1-D\nG1-G8"), (p, "2-D FD\nH1-H6")]):
        v = frame[col].to_numpy(float)
        jitter = (np.arange(v.size) - (v.size - 1) / 2) * 0.05
        ax.plot(i + jitter, v, "o", color=LINE, ms=5, mec=SURFACE, mew=1)
        for x, y, name in zip(i + jitter, v, frame.model):
            ax.annotate(name, (x, y), xytext=(4, 0), textcoords="offset points", fontsize=6.5, color=MUTED, va="center")
    ax.set_xticks([0, 1], ["Kennett 1-D\nG1-G8", "2-D FD\nH1-H6"])
    ax.set_xlim(-0.6, 1.8)
    ax.set_ylabel(ylab)
    tidy(ax, title)
axes[3].axhline(0, color=MUTED, lw=0.8)
axes[4].axhline(1, color=MUTED, lw=0.8)
axes[4].axhspan(1 / 1.1145, 1.1145, color=GRID, alpha=0.6, lw=0)
fig.text(0.01, 0.01, "(a-c) Standard MSE (R = 6) over the 14 factors, divided by that model's minimum; curves above 2.0 are clipped. "
         "(e) shaded: within one QScan grid step of Q0.", color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fd_curves.{ext}", dpi=200, facecolor=SURFACE)
print("written", FIG / "fd_curves.png")

# Gathers of the model with the largest throw.
name = "H4"   # 40 m throw
T = loadmat(ROOT / "truth" / f"{name}_truth.mat")["target"]
D = loadmat(ROOT / "models" / f"{name}.mat")["clean"]
t = np.arange(T.shape[0]) * 0.004
div = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e34948", "#8f1d1d"])
clip = np.percentile(np.abs(T), 99.5)
fig, axes = plt.subplots(1, 3, figsize=(11, 4.6), sharey=True, facecolor=SURFACE)
for ax, X, title in zip(axes, [T, D, D - T], ["Elastic target", "Noise-free attenuated data", "Attenuated - elastic"]):
    ax.imshow(X, aspect="auto", cmap=div, vmin=-clip, vmax=clip, extent=(0.5, X.shape[1] + 0.5, t[-1], t[0]), interpolation="none")
    ax.set_title(f"{name}: {title}", color=INK, fontsize=9, loc="left")
    ax.set_xlabel("trace")
axes[0].set_ylabel("two-way time (s)")
fig.text(0.01, 0.01, "Common colour scale (99.5th percentile of |target|); 2-D viscoacoustic FD, zero offset, 12 m traces.", color=MUTED, fontsize=7.5)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(FIG / "fd_gathers.png", dpi=200, facecolor=SURFACE)
print("written", FIG / "fd_gathers.png")
