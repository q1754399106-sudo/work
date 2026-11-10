"""Tables of pilot_mw.py against pilot_fd.py (same noise seeds; model means over realisations; descriptive).
usage: analyze_pilot_mw.py LEVEL"""
import sys, json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; level = sys.argv[1]


def ld(sub):
    return {(d["model"], d["r"]): d for d in (json.load(open(p, encoding="utf-8"))
            for p in glob.glob(str(ROOT / "results" / sub / f"*_{level}_r*.json")))}


M, F = ld("pilot_mw"), ld("pilot_fd"); keys = sorted(set(M) & set(F)); models = sorted({k[0] for k in keys})
ms = [r["m"] for r in M[keys[0]]["records"]]; i1 = ms.index(1.0)
KS = [k for k in M[keys[0]]["records"][0] if k not in ("m", "q_med")]; SC = sorted(k for k in KS if k.startswith("SC"))
tab = lambda src, k: np.array([[np.mean([src[q]["records"][j][k] for q in keys if q[0] == h]) for j in range(len(ms))] for h in models])
T = {k: tab(M, k) for k in KS}; T.update({k: tab(F, k) for k in ("UQW", "W0", "A0")})
ora = np.array([np.mean([F[q]["ORA"] for q in keys if q[0] == h]) for h in models])
print(f"{level}: {len(keys)} datasets; single-trace ceiling ORA per model", np.round(ora, 3), "mean %.3f" % ora.mean())
print("\nmean over models, by m:\n method    " + "  ".join("%5.2f" % m for m in ms) + "    mean")
for k in ["W0", "UQW", "A0"] + KS:
    v = T[k].mean(axis=0); print(" %-8s  " % k + "  ".join("%.3f" % x for x in v) + "   %.3f" % v.mean())
print("\nper model, scenario mean:\n method    " + "  ".join("%5s" % h for h in models))
for k in ["UQW"] + KS:
    print(" %-8s  " % k + "  ".join("%.3f" % x for x in T[k].mean(axis=1)))
print(" ORA(m=1)  " + "  ".join("%.3f" % x for x in ora))
pairs = [("UQMW_K9", "UQW"), ("MW0_K9", "UQW"), ("UQMW_K5", "UQW"), ("UQMW_K9", "MW0_K9"), ("UQMW_K1", "UQW")] + [("UQMW_K9", s) for s in SC]
for a, b in pairs:
    d = T[a].mean(axis=1) - T[b].mean(axis=1)
    print("scenario mean %s - %s: %+.3f  positive %d/%d  min %+.3f" % (a, b, d.mean(), (d > 0).sum(), len(d), d.min()))
d = T["MW0_K9"][:, i1] - ora; print("m = 1: MW0_K9 - ORA: %+.3f  positive %d/%d" % (d.mean(), (d > 0).sum(), len(d)))
d = T["UQMW_K9"].mean(axis=1) - ora; print("UQMW_K9 scenario mean - ORA(m=1): %+.3f  positive %d/%d" % (d.mean(), (d > 0).sum(), len(d)))
