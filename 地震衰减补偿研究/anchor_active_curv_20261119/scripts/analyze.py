"""Scores of PLAN_ZH.md: alpha_a = a0 sum(phi a) / sum(phi a^2) with phi = 0.1 s moving average of e(t) (judged) or
n(t) (reproduction of probe_active, not judged); per-profile median over datasets, errors against the published K.
-> analysis/active_curv.json, printed table.
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]; RES = HERE / "results"; OUT = HERE / "analysis"
WIN = 25                                                     # samples (4 ms) = 0.1 s, fixed in PLAN_ZH.md
# profile: (r, uniform alpha, K_ref) as in PLAN_ZH.md section 2
REF = {"base": (0.60, 0.675, 0.664), "S": (1.00, 0.880, 0.856), "D": (4.65, 0.294, 0.598), "D2": (2.28, 0.489, 0.624),
       "M": (0.24, 0.595, 0.614), "R1": (1.00, 0.628, 0.649), "R2": (0.42, 0.696, 0.689), "R3": (0.28, 0.598, 0.608),
       "R4": (0.38, 0.772, 0.737), "X1": (1.64, 0.556, 0.638), "X2": (3.96, 0.342, 0.611), "X3": (1.39, 0.619, 0.615)}
TRUNC = {2.2: (0.749, 0.80), 2.6: (0.448, 0.74)}           # D cut windows: (uniform alpha, K) from analysis/dtrunc.json
R_LE_1 = ["base", "S", "M", "R1", "R2", "R3", "R4"]


def smooth(x, w=WIN):
    x = np.asarray(x, float); k = np.ones(w)
    return np.convolve(x, k, "same") / np.convolve(np.ones_like(x), k, "same")


def alpha(rec, key):
    phi = smooth(rec[key]); a = np.asarray(rec["a"], float)
    return float(rec["a0"] * np.sum(phi * a) / np.sum(phi * a ** 2))


def uniform(rec):
    a = np.asarray(rec["a"], float); return float(rec["a0"] * a.mean() / (a ** 2).mean())


def main():
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(RES.glob("*.json"))]
    groups = {}
    for r in recs:
        groups.setdefault((r["profile"], r["t1"]), []).append(r)
    out = dict(window_samples=WIN, profiles={}, trunc={}, verdicts={})
    for (p, t1), rs in sorted(groups.items(), key=lambda kv: (kv[0][1] or 0, kv[0][0])):
        aE = [alpha(r, "e") for r in rs]; aN = [alpha(r, "n") for r in rs]; aU = [uniform(r) for r in rs]
        row = dict(n=len(rs), n_not_converged=sum(not r["converged"] for r in rs),
                   alpha_E_median=float(np.median(aE)), alpha_E_range=[min(aE), max(aE)],
                   alpha_N_median=float(np.median(aN)), alpha_N_range=[min(aN), max(aN)],
                   alpha_U_check=float(np.median(aU)), per_dataset={r["id"]: dict(E=e, N=nn) for r, e, nn in zip(rs, aE, aN)})
        if t1 is None:
            rr, aU0, K = REF[p]
            row.update(r=rr, alpha_U=aU0, K_ref=K, err_U=aU0 - K, err_E=row["alpha_E_median"] - K,
                       err_N=row["alpha_N_median"] - K, valid=row["n_not_converged"] <= len(rs) / 4)
            out["profiles"][p] = row
        else:
            aU0, K = TRUNC[t1]
            row.update(alpha_U=aU0, K_ref=K, err_U=aU0 - K, err_E=row["alpha_E_median"] - K, err_N=row["alpha_N_median"] - K)
            out["trunc"][str(t1)] = row
    P = out["profiles"]
    a_rows = {p: dict(abs_err_E=abs(P[p]["err_E"]), limit=abs(P[p]["err_U"]) + 0.02,
                      ok=bool(P[p]["valid"] and abs(P[p]["err_E"]) <= abs(P[p]["err_U"]) + 0.02)) for p in R_LE_1}
    out["verdicts"]["a"] = dict(rows=a_rows, holds=all(v["ok"] for v in a_rows.values()))
    eD = P["D"]["err_E"]
    out["verdicts"]["b"] = dict(abs_err_E_D=abs(eD), limit=0.5 * 0.304, holds=bool(P["D"]["valid"] and abs(eD) <= 0.152))
    OUT.mkdir(exist_ok=True); (OUT / "active_curv.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"{'profile':8s} {'r':>5s} {'K':>6s} {'U':>6s} {'eU':>7s} | {'E med':>6s} [range]       {'eE':>7s} | {'N med':>6s} {'eN':>7s}  nc")
    for p, v in P.items():
        print(f"{p:8s} {v['r']:5.2f} {v['K_ref']:6.3f} {v['alpha_U']:6.3f} {v['err_U']:+7.3f} | {v['alpha_E_median']:6.3f} "
              f"[{v['alpha_E_range'][0]:.3f},{v['alpha_E_range'][1]:.3f}] {v['err_E']:+7.3f} | {v['alpha_N_median']:6.3f} "
              f"{v['err_N']:+7.3f}  {v['n_not_converged']}  (U check {v['alpha_U_check']:.3f})")
    for t1, v in out["trunc"].items():
        print(f"D t1={t1}  K {v['K_ref']:.3f} U {v['alpha_U']:.3f} eU {v['err_U']:+.3f} | E {v['alpha_E_median']:.3f} "
              f"eE {v['err_E']:+.3f} | N {v['alpha_N_median']:.3f} eN {v['err_N']:+.3f}  nc {v['n_not_converged']}")
    print("a:", out["verdicts"]["a"]["holds"], {p: v["ok"] for p, v in a_rows.items()})
    print("b:", out["verdicts"]["b"])


if __name__ == "__main__":
    main()
