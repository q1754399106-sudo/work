"""Paper figure 1 -> figures/schematic.png: cumulative t* and kappa against the wavelet-window centre.

(a) a(t) of the three profiles in the analysis window, the wavelet-window centre a0 (1.4 s) and the weighted centre
<a^2>/<a> of the base profile. (b) kappa = a0<a>/<a^2> as a function of the wavelet-window centre for the three
profiles (uniform weights), with the measured anchor coefficients (pre-registered tests; D in red).
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(Path(__file__).parent))
from env15 import lib, semi  # noqa: E402

plt.rcParams.update({"font.sans-serif": ["Noto Sans CJK SC", "WenQuanYi Zen Hei", "DejaVu Sans"],
                     "axes.unicode_minus": False, "font.size": 9})
cm = semi.curtain_cm(); t = semi.TW; tt = np.linspace(0, 3.0, 601)
PROF = {"基准 78.5/78.5/135/135": ([78.5, 78.5, 135, 135], "k"), "S 30/200/200/200": ([30, 200, 200, 200], "#2a6fbb"),
        "D 250/250/40/40": ([250, 250, 40, 40], "#b8322e")}
an = json.loads((ROOT / "analysis" / "anchor.json").read_text())["by_weight"]
win = json.loads((ROOT / "analysis" / "window.json").read_text())["windows"]
pr = json.loads((ROOT / "analysis" / "profile.json").read_text())

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for name, (Q, col) in PROF.items():
    c = dict(cm, Q=np.asarray(Q, float)); a = lib.absorption(c, tt); ax[0].plot(tt, a, color=col, label=name)
    aw = lib.absorption(c, t); tc = np.linspace(1.25, 2.95, 60)
    k = [lib.absorption(c, [x])[0] * aw.mean() / (aw ** 2).mean() for x in tc]; ax[1].plot(tc, k, color=col, label=name)
cb = dict(cm, Q=np.asarray(PROF["基准 78.5/78.5/135/135"][0], float)); ab = lib.absorption(cb, t)
ax[0].axvspan(1.2, 3.0, color="#f2f2f2", zorder=0); ax[0].axvspan(1.2, 1.6, color="#dde8f5", zorder=0)
a0 = lib.absorption(cb, [1.4])[0]; ax[0].axhline(a0, color="k", ls=":", lw=0.8)
ax[0].axhline((ab ** 2).mean() / ab.mean(), color="k", ls="--", lw=0.8)
ax[0].text(0.05, a0, " a₀（子波窗 1.2–1.6 s 中点）", va="bottom", fontsize=7)
ax[0].text(0.05, (ab ** 2).mean() / ab.mean(), " ⟨a²⟩/⟨a⟩（分析窗加权中心）", va="bottom", fontsize=7)
ax[0].set_xlabel("双程时 t (s)"); ax[0].set_ylabel("累积 t*  a(t) (s)"); ax[0].legend(fontsize=7, loc="upper left")
ax[0].set_title("(a) 吸收剖面与窗口（浅灰：分析窗，浅蓝：子波窗）")
ax[1].axhline(1, color="#999", lw=0.8); ax[1].axhline(0, color="#999", lw=0.8)
meas = [(1.4, an["0.5"]["K"], "k"), (1.4, an["0.0625"]["K"], "k"), (2.0, win["2.0"]["K"], "k"), (2.6, win["2.6"]["K"], "k"),
        (1.4, pr["S"]["K"], "#2a6fbb"), (1.4, pr["D"]["K"], "#b8322e")]
for x, y, col in meas:
    ax[1].plot(x, y, "o", color=col, ms=6, mec="white")
ax[1].text(1.45, pr["D"]["K"] + 0.03, "D 实测（预测失败）", color="#b8322e", fontsize=7)
ax[1].text(1.27, 1.03, "κ > 1：选择与真 Q 反向", fontsize=7, color="#666")
ax[1].set_xlabel("子波窗中点 (s)"); ax[1].set_ylabel("κ = a₀⟨a⟩/⟨a²⟩"); ax[1].legend(fontsize=7, loc="lower right")
ax[1].set_title("(b) κ 随子波窗深度与剖面形状（点：实测）")
fig.tight_layout(); fig.savefig(ROOT / "figures" / "schematic.png", dpi=160); print("figures/schematic.png")
