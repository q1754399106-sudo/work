"""Score PLAN_CRIT_ZH.md from results17/ -> analysis/crit.json."""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; DV = 1 / 0.7 - 1 / 1.4
R = {}
for f in glob.glob(str(ROOT / "results17" / "*.json")):
    r = json.loads(Path(f).read_text()); R.setdefault(r["dataset"], {})[r["c_anchor"]] = r
out = {}
for crit in ("J", "chengD", "yan", "lp1", "lp05"):
    o = {}
    for kind in ("S", "F"):
        ds = {d: v for d, v in R.items() if d.startswith(kind + "_") and 0.7 in v and 1.4 in v}
        recs = [r for v in ds.values() for r in v.values()]
        edge = sum(r["sel"][crit]["edge"] for r in recs) / max(1, len(recs))
        bad = sum(not r["sel"][crit]["ok"] for r in recs) / max(1, len(recs))
        k = [(v[0.7]["sel"][crit]["v_sel"] - v[1.4]["sel"][crit]["v_sel"]) / DV for v in ds.values()]
        e = dict(n=len(ds), edge_share=edge, unconverged_share=bad, kappa_median=float(np.median(k)) if k else None,
                 kappa_range=[float(min(k)), float(max(k))] if k else None, kappa=dict(zip(ds, map(float, k))))
        if kind == "S":
            rows = [r for v in ds.values() for r in v.values()]
            X = np.array([[r["u"], r["v_anchor"], 1.0] for r in rows]); y = np.array([r["sel"][crit]["v_sel"] for r in rows])
            A, K, C = np.linalg.lstsq(X, y, rcond=None)[0]; e.update(A=float(A), K=float(K), C=float(C), A_plus_K=float(A + K))
        o[kind] = e
    out[crit] = o
V = {}
ok = lambda c, k: out[c][k]["unconverged_share"] <= 0.25
V["K-J"] = {k: bool(ok("J", k) and abs(out["J"][k]["kappa_median"] - 0.675) <= 0.12) for k in ("S", "F")}
kS, kF = out["chengD"]["S"]["kappa_median"], out["chengD"]["F"]["kappa_median"]
if ok("chengD", "S") and ok("chengD", "F"):
    V["K-C"] = "anchored" if (kS > 0.4 and kF > 0.4) else ("immune" if (kS < 0.2 and kF < 0.2) else "undecided")
else:
    V["K-C"] = "undecided (convergence)"
V["K-Y"] = {k: ("degenerate" if out["yan"][k]["edge_share"] > 0.25 else
               ("anchored" if out["yan"][k]["kappa_median"] > 0.4 else "not anchored")) for k in ("S", "F")}
V["lp_edge_share"] = {c: {k: out[c][k]["edge_share"] for k in ("S", "F")} for c in ("lp1", "lp05")}
out["verdicts"] = V
(ROOT / "analysis" / "crit.json").write_text(json.dumps(out, indent=1))
for c in ("J", "chengD", "yan", "lp1", "lp05"):
    for k in ("S", "F"):
        print(c, k, {x: (round(y, 3) if isinstance(y, float) else y) for x, y in out[c][k].items() if x != "kappa"})
print(V)
