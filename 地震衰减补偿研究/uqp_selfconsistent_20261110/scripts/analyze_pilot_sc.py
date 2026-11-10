"""SC baseline (pilot_sc.py) against UQ-W (pilot_fd.py) and UQ-MW (pilot_mw.py) on the same noisy gathers.
Descriptive. usage: analyze_pilot_sc.py LEVEL"""
import sys, json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; level = sys.argv[1]
ld = lambda sub: {(d["model"], d["r"]): d for d in (json.load(open(p, encoding="utf-8"))
                  for p in glob.glob(str(ROOT / "results" / sub / f"*_{level}_r*.json")))}
S, M, F = ld("pilot_sc"), ld("pilot_mw"), ld("pilot_fd"); keys = sorted(set(S) & set(M) & set(F)); ms = [r["m"] for r in S[keys[0]]["records"]]
g = lambda src, q, m, k: [r for r in src[q]["records"] if r["m"] == m][0][k]
rows = [("SC_orig", S), ("SC_tuned", S), ("UQW", F), ("MW0_K9", M), ("UQMW_K9", M)]
T = {k: np.array([[g(src, q, m, k) for m in ms] for q in keys]) for k, src in rows}
print(f"{level}: datasets", ["%s r%d" % q for q in keys])
print("mean over datasets, by m:\n method     " + "  ".join("%5.2f" % m for m in ms) + "    mean")
for k, _ in rows:
    v = T[k].mean(axis=0); print(" %-9s  " % k + "  ".join("%.3f" % x for x in v) + "   %.3f" % v.mean())
print("per dataset (mean over m):\n method     " + "  ".join("%5s" % q[0] for q in keys))
for k, _ in rows:
    print(" %-9s  " % k + "  ".join("%.3f" % x for x in T[k].mean(axis=1)))
for a, b in (("UQMW_K9", "SC_orig"), ("MW0_K9", "SC_orig"), ("SC_orig", "UQW")):
    d = T[a].mean(axis=1) - T[b].mean(axis=1); print("%s - %s: %+.3f  positive %d/%d  min %+.3f" % (a, b, d.mean(), (d > 0).sum(), len(d), d.min()))
i1 = ms.index(1.0); d = T["MW0_K9"][:, i1] - T["SC_orig"][:, i1]; print("m = 1: MW0_K9 - SC_orig: %+.3f  positive %d/%d" % (d.mean(), (d > 0).sum(), len(d)))
conv = lambda k: sum(r[k + "_conv"] for q in keys for r in S[q]["records"]); n = len(keys) * len(ms)
print("converged: SC_orig %d/%d, SC_tuned %d/%d;  median seconds per solve: SC_orig %.0f, SC_tuned %.0f" % (
    conv("SC_orig"), n, conv("SC_tuned"), n, np.median([r["SC_orig_s"] for q in keys for r in S[q]["records"]]),
    np.median([r["SC_tuned_s"] for q in keys for r in S[q]["records"]])))
