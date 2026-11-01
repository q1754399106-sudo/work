"""Field figure -> figures/field.png: two-anchor kappa, implied beta, and the current QScan vs SR along the survey."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.sans-serif": ["Noto Sans CJK SC", "WenQuanYi Zen Hei", "DejaVu Sans"],
                     "axes.unicode_minus": False, "font.size": 9})
F2 = json.loads((ROOT / "analysis" / "field2.json").read_text()); F1 = json.loads((ROOT / "analysis" / "field.json").read_text())
L = sorted(int(k) for k in F2["lines"]); g = lambda k: np.array([F2["lines"][str(x)][k] for x in L])
L1 = sorted(int(k) for k in F1["lines"])

fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))
ax[0].plot(L, g("kappa"), "o-", color="#2e8b57", ms=4, label="26 条 SR 测线")
ax[0].plot(L1, [F1["lines"][str(x)]["kappa"] for x in L1], "s", mfc="none", color="#2e8b57", label="第 5j 节 10 条测线")
ax[0].axhline(0.675, color="k", ls="--", lw=0.8, label="模型 Q 算出的 κ = 0.675"); ax[0].set_ylim(0.5, 0.85)
ax[0].set_xlabel("inline"); ax[0].set_ylabel("两锚点 κ̂"); ax[0].legend(fontsize=7); ax[0].set_title("(a) κ̂ 在全工区不变")
ax[1].plot(L, g("beta"), "o-", color="#b8322e", ms=4); ax[1].axhspan(0.03, 0.07, color="#ddd", label="半合成上的 β")
ax[1].set_ylim(0, 0.4); ax[1].set_xlabel("inline"); ax[1].set_ylabel("β̂（以 SR 为准）"); ax[1].legend(fontsize=7)
ax[1].set_title("(b) 实测截距 β 约 0.28，比半合成大")
cs = g("c_sr"); c1 = g("c1"); ax[2].plot(cs, c1, "o", color="#2a6fbb", ms=4, label="现流程 QScan（顶点）")
xx = np.linspace(0.8, 1.4, 2); k, b = np.median(g("kappa")), np.median(g("beta"))
ax[2].plot(xx, 1 / ((1 - k) / xx + k + b), color="#2a6fbb", lw=1, ls="--", label="1/[(1−κ̂)/c + κ̂ + β̂]")
ax[2].plot(xx, xx, color="k", lw=0.8, label="c_QScan = c_SR")
ax[2].set_xlabel("c_SR"); ax[2].set_ylabel("c_QScan"); ax[2].legend(fontsize=7); ax[2].set_title("(c) QScan 几乎不随 SR 变")
fig.tight_layout(); fig.savefig(ROOT / "figures" / "field.png", dpi=150); print("figures/field.png")
