"""Summarise the beta ladder (BETA_ZH.md) from results12/ -> analysis/beta_ladder.json."""
import json, glob, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; DV = 1 / 0.7 - 1 / 1.4; RNG = np.random.default_rng(20261107)
R = {}
for f in glob.glob(str(ROOT / "results12" / "*.json")):
    r = json.loads(Path(f).read_text()); R.setdefault(r["variant"], {}).setdefault((r["c_true"], r["seed"]), {})[r["c_anchor"]] = r
out = {}
for v, d in sorted(R.items()):
    rows = []
    for (c, s), x in sorted(d.items()):
        if len(x) < 3:
            continue
        k = (x[0.7]["v_sel"] - x[1.4]["v_sel"]) / DV; u = 1 / c
        rows.append(dict(c=c, seed=s, kappa=k, beta=x[1.0]["v_sel"] - (1 - k) * u - k, v1=x[1.0]["v_sel"],
                         bad=sum(r["edge"] or not r["neighbourhood_converged"] for r in x.values())))
    b = np.array([r["beta"] for r in rows]); bs = [b[RNG.integers(0, len(b), len(b))].mean() for _ in range(2000)]
    out[v] = dict(n=len(rows), beta_mean=float(b.mean()), beta_ci=[float(np.percentile(bs, 5)), float(np.percentile(bs, 95))],
                  kappa_mean=float(np.mean([r["kappa"] for r in rows])), unreliable=int(sum(r["bad"] for r in rows)),
                  by_c={str(c): float(np.mean([r["beta"] for r in rows if r["c"] == c])) for c in sorted({r["c"] for r in rows})},
                  rows=rows)
(ROOT / "analysis" / "beta_ladder.json").write_text(json.dumps(out, indent=1))
for v, s in out.items():
    print(v, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in s.items() if k != "rows"})
