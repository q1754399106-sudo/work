"""figures/mu_effects.png: per-dataset S differences of the mu comparisons (left) and mean S per method (right)."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
s = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8")); P = s["per_dataset"]
PAIRS = [("SR_mu0", "SR_mu10", "SR: mu 0 - 10 (H1)"), ("QT_mu0", "QT_mu10", "QScan, a=2^-13: mu 0 - 10 (H2)"),
         ("OR_mu0", "OR_mu10", "true Q: mu 0 - 10 (H3)"), ("SR_mu1", "SR_mu10", "SR: mu 1 - 10 (descriptive)"),
         ("SR_mu1", "SR_mu0", "SR: mu 1 - 0 (descriptive)")]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw=dict(width_ratios=[1.5, 1]))
cols = {"D": "tab:blue", "E": "tab:orange", "F": "tab:green"}
for i, (a, b, lab) in enumerate(PAIRS):
    d = np.array([p["S"][a] - p["S"][b] for p in P])
    for j, p in enumerate(P):
        a1.scatter(i + (j - 5.5) * 0.04, d[j], color=cols[p["id"][1]], marker="o" if p["id"].find("0.75") > 0 else "^", s=22)
    a1.hlines(d.mean(), i - 0.3, i + 0.3, color="k")
a1.axhline(0, color="grey", lw=0.6); a1.set_xticks(range(len(PAIRS))); a1.set_xticklabels([p[2] for p in PAIRS], rotation=20, ha="right", fontsize=8)
a1.set_ylabel("paired difference of S"); a1.set_title("12 sealed datasets (colour: structure D/E/F; o c0 0.75, ^ 1.05)", fontsize=9)
M = ["IQ70", "VDB", "SR_mu10", "SR_mu1", "SR_mu0", "QSCa", "QT_mu10", "QT_mu0", "OR_mu10", "OR_mu0"]
a2.barh(range(len(M)), [s["mean_S"][m] for m in M], color=["0.6", "0.6", "tab:red", "tab:red", "tab:red", "tab:purple",
                                                          "tab:purple", "tab:purple", "tab:cyan", "tab:cyan"])
for k, m in enumerate(M):
    a2.text(s["mean_S"][m] + 0.005, k, f"{s['mean_S'][m]:.3f}", va="center", fontsize=8)
a2.set_yticks(range(len(M))); a2.set_yticklabels(M, fontsize=8); a2.set_xlim(0.6, 0.95); a2.invert_yaxis(); a2.set_xlabel("mean S")
fig.tight_layout(); (ROOT / "figures").mkdir(exist_ok=True); fig.savefig(ROOT / "figures" / "mu_effects.png", dpi=110)
