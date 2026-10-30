"""figures/scores.png (mean S per method, Z and Y; local windows of Z) and figures/compression.png (QScan 1/c
against the theory-weighted u_QS, Z and Y) from results/test and analysis/summary.json. Descriptive."""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8"))
M = ["Input", "IQ30", "IQ70", "SCa_SR", "VDB", "QSCa", "QSCa_tuned", "SCa_oracle", "IQ70_oracle"]
FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)

fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
x = np.arange(len(M)); w = 0.38
axs[0].bar(x - w / 2, [S["Z"]["mean_S"][m] for m in M], w, label="Z (complex Q)", color="#3b6ea5")
axs[0].bar(x + w / 2, [S["Y"]["mean_S"][m] for m in M], w, label="Y (twin, proportional Q)", color="#b9c9de")
axs[0].set_xticks(x, M, rotation=40, ha="right"); axs[0].set_ylabel("mean S (B, C; 12 each)"); axs[0].legend(fontsize=8)
axs[0].set_ylim(0, 1)
K = ["far", "fault", "beside_gas", "below_gas"]
for i, k in enumerate(K):
    axs[1].bar(x + (i - 1.5) * 0.2, [S["Z"]["local_mean"][m][k] for m in M], 0.2, label=k)
axs[1].set_xticks(x, M, rotation=40, ha="right"); axs[1].set_ylabel("pooled correlation, Z, 2.0-3.0 s")
axs[1].legend(fontsize=8); axs[1].set_ylim(-0.4, 1)
fig.tight_layout(); fig.savefig(FIG / "scores.png", dpi=120); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.5, 5))
for g, col in (("Z", "#3b6ea5"), ("Y", "#c0504d")):
    sel = S[g]["selections"]; u = [s["u_qs"] for s in sel]; v = [1 / s["qscan"] for s in sel]
    ax.scatter(u, v, color=col, label=f"{g}: A = {S[g]['H4']['A']:.2f} [{S[g]['H4']['lo']:.2f}, {S[g]['H4']['hi']:.2f}]")
uu = np.linspace(0.6, 1.6, 2)
ax.plot(uu, uu, "k:", lw=1, label="slope 1")
ax.plot(uu, 0.325 * uu + 0.675, "k--", lw=1, label="theory 1 - kappa = 0.33 (beta = 0)")
ax.set_xlabel("u_QS (theory-weighted true 1/c)"); ax.set_ylabel("1 / c selected by QScan"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "compression.png", dpi=120); plt.close(fig)
