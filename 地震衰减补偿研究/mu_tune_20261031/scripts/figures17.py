"""figures/mu_tune.png: development mean S over the mu grid (left) and the confirmation paired differences (right)."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "field_application_20261028" / "cache" / "deps"))   # matplotlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
dev = json.loads((ROOT / "analysis" / "dev_summary.json").read_text(encoding="utf-8"))
s = json.loads((ROOT / "analysis" / "summary.json").read_text(encoding="utf-8")); P = s["per_dataset"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8))
x = np.arange(7); labels = list(dev["mean_S"]["SR"])
for route, col, lab in (("SR", "tab:red", "SR factor"), ("QT", "tab:purple", "QScan factor")):
    a1.plot(x, [dev["mean_S"][route][m] for m in labels], "o-", color=col, label=f"SC-a (a = 2^-13), {lab}")
a1.axhline(dev["IQ70"], color="0.5", ls="--", label="IQ70 (SR factor)")
a1.set_xticks(x); a1.set_xticklabels(labels); a1.set_xlabel("mu_relative"); a1.set_ylabel("mean S")
a1.set_title("development: 24 used datasets (A-F)", fontsize=9); a1.legend(fontsize=8)
PAIRS = [("SR_sel", "SR_mu10", "C1 SR: 0.25 - 10"), ("QT_sel", "QT_mu10", "C2 QScan: 0.25 - 10"),
         ("SR_sel", "IQ70", "C3 SR 0.25 - IQ70"), ("QT_sel", "QSCa", "C4 user total"), ("QT_mu10", "QSCa", "only smaller a")]
cols = {"G": "tab:blue", "H": "tab:orange", "I": "tab:green"}
for i, (a, b, lab) in enumerate(PAIRS):
    d = np.array([p["S"][a] - p["S"][b] for p in P])
    for j, p in enumerate(P):
        a2.scatter(i + (j - 5.5) * 0.04, d[j], color=cols[p["id"][1]], marker="o" if "0.75" in p["id"] else "^", s=22)
    a2.hlines(d.mean(), i - 0.3, i + 0.3, color="k")
a2.axhline(0, color="grey", lw=0.6); a2.set_xticks(range(len(PAIRS))); a2.set_xticklabels([p[2] for p in PAIRS], rotation=15, fontsize=8)
a2.set_ylabel("paired difference of S"); a2.set_title("confirmation: 12 sealed datasets (colour G/H/I; o c0 0.75, ^ 1.05)", fontsize=9)
fig.tight_layout(); (ROOT / "figures").mkdir(exist_ok=True); fig.savefig(ROOT / "figures" / "mu_tune.png", dpi=110)
