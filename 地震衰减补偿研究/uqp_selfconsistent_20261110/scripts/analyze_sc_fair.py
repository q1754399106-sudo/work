"""SC-a fair-control tables at 4 x noise (pilot_sc_tune.py eval, pilot_sc_robust.py) against SR -> MW0 K 9 and W0
(pilot_semi.py). 16 evaluation datasets; 90% bootstrap intervals, 2000 draws. Descriptive."""
import json, glob, collections
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; rng = np.random.default_rng(20261125); rd = lambda p: json.loads(Path(p).read_text(encoding="utf-8"))


def ci(d):
    b = rng.integers(0, len(d), (2000, len(d))); m = d[b].mean(axis=1)
    return "%+.3f [%+.3f, %+.3f] pos %d/%d" % (d.mean(), np.quantile(m, .05), np.quantile(m, .95), (d > 0).sum(), len(d))


V = collections.defaultdict(dict)
for p in glob.glob(str(ROOT / "results/pilot_sc_robust/*.json")):
    r = rd(p); V[r["variant"]][r["id"]] = r
for p in glob.glob(str(ROOT / "results/pilot_sc_tune/eval_*.json")):
    r = rd(p); V["R0_sel"][r["id"]] = r
R = {d["id"]: d for d in map(rd, glob.glob(str(ROOT / "results/pilot_semi/*_nd.json")))}
for v, D in sorted(V.items()):
    ids = sorted(D); sc = np.array([D[i]["S"] for i in ids]); mw = np.array([R[i]["sr"]["MW0_K9"] for i in ids]); w0 = np.array([R[i]["sr"]["W0"] for i in ids])
    fs = np.median([D[i].get("fmax_shallow", np.nan) for i in ids]); fd = np.median([D[i].get("fmax_target", np.nan) for i in ids])
    print("%-7s n %d  a 2^%+.0f mu %g  SC %.3f  MW0_K9 %.3f  W0 %.3f  fs %.0f fd %.0f  conv %d  median s %.0f" % (
        v, len(ids), np.log2(D[ids[0]]["a"]), D[ids[0]]["mu"], sc.mean(), mw.mean(), w0.mean(), fs, fd, sum(D[i]["converged"] for i in ids), np.median([D[i]["seconds"] for i in ids])))
    print("        MW0_K9 - SC", ci(mw - sc), "   SC - W0", ci(sc - w0))
