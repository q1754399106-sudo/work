"""Post hoc (not in PLAN_ITER_ZH.md): increment ratio of the mid / deep arms using only reliable steps.

A ratio Delta_k / Delta_{k-1} is used only when steps k-1 and k both had a converged neighbourhood (selection and
its two grid neighbours converged within 20000 FISTA iterations) and neither hit the grid edge. Pooled
least-squares ratio as in analyze_iter.py, bootstrap over datasets (2000, seed 20261103).
-> analysis/iter_posthoc.json
"""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RNG = np.random.default_rng(20261103)
KAPPA = dict(main=0.675, mid=0.964, deep=1.140)


def pairs(r):
    h = r["history"]; v = [1.0] + [x["v_sel"] for x in h]; ok = [x["neighbourhood_converged"] and not x["edge"] for x in h]
    d = np.diff(v); return [(d[k - 1], d[k]) for k in range(1, len(d)) if ok[k - 1] and ok[k]]


def ratio(recs):
    P = [p for r in recs for p in pairs(r)]
    if not P:
        return float("nan"), 0
    a = np.array(P); return float(np.sum(a[:, 0] * a[:, 1]) / np.sum(a[:, 0] ** 2)), len(P)


out = {}
for arm in ("main", "mid", "deep"):
    recs = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(ROOT / "results5" / f"*__{arm}.json")))]
    if not recs:
        continue
    r, n = ratio(recs); bs = [ratio([recs[i] for i in RNG.integers(0, len(recs), len(recs))])[0] for _ in range(2000)]
    out[arm] = dict(kappa=KAPPA[arm], n_datasets=len(recs), n_pairs=n, r_reliable=r,
                    r_ci=[float(np.nanpercentile(bs, 5)), float(np.nanpercentile(bs, 95))],
                    per_dataset={x["id"]: [round(b / a, 3) for a, b in pairs(x)] for x in recs})
(ROOT / "analysis" / "iter_posthoc.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for arm, s in out.items():
    print(arm, {k: v for k, v in s.items() if k != "per_dataset"})
