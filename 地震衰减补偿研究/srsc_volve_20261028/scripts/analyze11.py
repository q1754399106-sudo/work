"""Verdicts of experiment 11 (PLAN_ZH.md section 4) -> analysis/summary.json."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sc11  # noqa: E402
from sc11 import np  # noqa: E402

T90 = {7: 1.895, 19: 1.729}
MAIN, BASE = "SRSCa_1.3+ph", "QSCa"


def ci(d):
    d = np.asarray(d, float); n = len(d); m = float(d.mean()); se = float(d.std(ddof=1) / np.sqrt(n))
    return dict(n=n, mean=m, ci90=[m - T90[n - 1] * se, m + T90[n - 1] * se], positive=int((d > 0).sum()))


R = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((sc11.ROOT / "results" / "synth").glob("*.json"))]
G = {g: [r for r in R if r["group"] == g] for g in ("T", "L")}
assert len(G["T"]) == 20 and len(G["L"]) == 8, {g: len(v) for g, v in G.items()}
methods = list(G["T"][0]["scores"]); out = dict(means={}, vs_QSCa={}, sr={})
for g, rs in G.items():
    out["means"][g] = {m: {k: float(np.mean([r["scores"][m][k] for r in rs])) for k in ("M1", "M1_lag")} for m in methods}
    out["vs_QSCa"][g] = {m: {k: ci([r["scores"][m][k] - r["scores"][BASE][k] for r in rs]) for k in ("M1", "M1_lag")}
                         for m in methods if m != BASE}
    e = np.array([np.log(r["sr"]["c"] / r["c_true"]) for r in rs]); q = np.array([np.log(r["qscan"] / r["c_true"]) for r in rs])
    out["sr"][g] = dict(mean_log=float(e.mean()), mean_abs_log=float(np.abs(e).mean()), within25=float(np.mean(np.abs(e) <= np.log(1.25))),
                        qscan_mean_abs_log=float(np.abs(q).mean()))
    out["means"][g + "_by_c"] = {str(c): {m: float(np.mean([r["scores"][m]["M1"] for r in rs if r["c_true"] == c])) for m in methods}
                                 for c in sorted({r["c_true"] for r in rs})}
t, l = out["vs_QSCa"]["T"][MAIN], out["vs_QSCa"]["L"][MAIN]
out["P1"] = bool(t["M1"]["mean"] > 0 and t["M1"]["ci90"][0] > 0)
out["P2"] = bool(t["M1_lag"]["ci90"][0] > -0.02)
out["P3"] = bool(l["M1"]["ci90"][0] > -0.05 and l["M1_lag"]["ci90"][0] > -0.05)
(sc11.ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print({k: out[k] for k in ("P1", "P2", "P3")})
for g in ("T", "L"):
    print(g, "SR", {k: round(v, 3) for k, v in out["sr"][g].items()})
    for m in methods:
        mm = out["means"][g][m]; d = out["vs_QSCa"][g].get(m)
        print(f"  {m:16s} M1 {mm['M1']:.3f} lag {mm['M1_lag']:.3f}",
              "" if d is None else f"| dM1 {d['M1']['mean']:+.3f} [{d['M1']['ci90'][0]:+.3f},{d['M1']['ci90'][1]:+.3f}] "
                                   f"dlag {d['M1_lag']['mean']:+.3f} [{d['M1_lag']['ci90'][0]:+.3f},{d['M1_lag']['ci90'][1]:+.3f}]")
