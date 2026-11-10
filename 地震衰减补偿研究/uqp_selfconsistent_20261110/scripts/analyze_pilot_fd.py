"""Tables of pilot_fd.py results (model means over realisations; descriptive). usage: analyze_pilot_fd.py LEVEL"""
import sys, json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; level = sys.argv[1]
KEYS = ("W0", "UQW", "A0", "UQA", "PW0", "PA0", "PUQW", "PUQ")
files = sorted(glob.glob(str(ROOT / "results" / "pilot_fd" / f"*_{level}_r*.json")))
R = [json.load(open(p, encoding="utf-8")) for p in files]; models = sorted({d["model"] for d in R})
ms = [r["m"] for r in R[0]["records"]]
tab = {k: np.array([[np.mean([d["records"][j][k] for d in R if d["model"] == h]) for j in range(len(ms))] for h in models]) for k in KEYS}
ora = np.array([np.mean([d["ORA"] for d in R if d["model"] == h]) for h in models])
orab = np.array([np.mean([max([d["ORA"]] + list(d["ORA_u"].values())) for d in R if d["model"] == h]) for h in models])
print(f"{level}: {len(R)} files, models {models}")
print("ceiling ORA (true Q) per model:", np.round(ora, 3), "mean %.3f;  best u of the grid: mean %.3f" % (ora.mean(), orab.mean()))
print("\nmean over models, by m:\n method  " + "  ".join("%5.2f" % m for m in ms) + "    mean   m=1-ORA")
for k in KEYS:
    v = tab[k].mean(axis=0); print(" %-6s  " % k + "  ".join("%.3f" % x for x in v) + "   %.3f   %+.3f" % (v.mean(), v[ms.index(1.0)] - ora.mean()))
print("\nper model, scenario mean (9 m):\n method  " + "  ".join("%5s" % h for h in models))
for k in KEYS:
    print(" %-6s  " % k + "  ".join("%.3f" % x for x in tab[k].mean(axis=1)))
print("\nper model at m = 1:\n method  " + "  ".join("%5s" % h for h in models))
for k in KEYS:
    print(" %-6s  " % k + "  ".join("%.3f" % x for x in tab[k][:, ms.index(1.0)]))
print(" ORA     " + "  ".join("%.3f" % x for x in ora))
uf = np.array([[np.mean([d["records"][j]["u_fit"] for d in R if d["model"] == h]) for j in range(len(ms))] for h in models])
se = np.array([[np.median([d["records"][j]["se"] for d in R if d["model"] == h]) for j in range(len(ms))] for h in models])
print("\nfitted u / m (truth 1) per model, mean over m:", np.round((uf / np.array(ms)[None, :]).mean(axis=1), 2), " jackknife se (median):", np.round(np.median(se, axis=1), 2))
for a, b in (("A0", "W0"), ("UQA", "UQW"), ("A0", "UQW"), ("PUQ", "UQW")):
    d = tab[a].mean(axis=1) - tab[b].mean(axis=1); print("scenario mean %s - %s: %+.3f  positive models %d/%d" % (a, b, d.mean(), (d > 0).sum(), len(d)))
