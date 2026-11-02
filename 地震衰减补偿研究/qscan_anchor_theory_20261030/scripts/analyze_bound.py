"""Score the validity-boundary test (PLAN_BOUND_ZH.md) from results26/ -> analysis/bound.json.
Per profile: two-anchor kappa (median over datasets), current-pipeline slope A_cur, A + K, bootstrap 90% intervals
(as analyze_profile2.py); verdicts B1-B3; descriptive Spearman correlation of r and (K - kappa_U)."""
import json, glob
from pathlib import Path
import numpy as np
import analyze_profile2 as ap2

ROOT = Path(__file__).resolve().parents[1]
# uniform kappa and deep/shallow rate ratio r = [a(3.0) - a(2.1)] / [a(2.1) - a(1.2)], computed before the run
KU = {"R1": 0.628, "R2": 0.6958, "R3": 0.5975, "R4": 0.7716, "X1": 0.5562, "X2": 0.3417, "X3": 0.619}
R = {"R1": 1.0, "R2": 0.416, "R3": 0.28, "R4": 0.376, "X1": 1.637, "X2": 3.96, "X3": 1.393}
OLD = {"base": (0.675, 0.60, 0.66), "S": (0.880, 1.00, 0.856), "D": (0.294, 4.65, 0.598), "D2": (0.489, 2.28, 0.624),
       "M": (0.595, 0.24, 0.614)}   # kappa_U, r, measured K (results7 / results10 / anchor experiment)


def spearman(x, y):
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y)); return float(np.corrcoef(rx, ry)[0, 1])


if __name__ == "__main__":
    recs = [json.loads(Path(f).read_text()) for f in glob.glob(str(ROOT / "results26" / "*.json"))]
    ap2.KU.clear(); ap2.KU.update(KU); ap2.RNG = np.random.default_rng(20261102)
    S = {p: ap2.summary(p, recs) for p in KU if any(r["profile"] == p for r in recs)}
    for p in S:
        S[p]["r"] = R[p]; S[p]["err"] = S[p]["kappa_median"] - KU[p]; S[p]["valid"] = S[p]["unreliable_share"] <= 0.25
    V = {}
    low = [p for p in S if R[p] <= 1.0 + 1e-6 and S[p]["valid"]]; high = [p for p in S if R[p] > 1.0 + 1e-6 and S[p]["valid"]]
    V["B1"] = dict(profiles=low, holds=bool(low) and all(abs(S[p]["err"]) <= 0.10 for p in low))
    V["B2"] = dict(profiles=high, holds=bool(high) and all(S[p]["err"] >= -0.05 for p in high))
    V["B3"] = dict(profiles=[p for p in S if S[p]["valid"]],
                   holds=all(abs(S[p]["A_current"] - S[p]["predicted_A"]) <= 0.15 for p in S if S[p]["valid"]))
    xs = [R[p] for p in S] + [v[1] for v in OLD.values()]; ys = [S[p]["err"] for p in S] + [v[2] - v[0] for v in OLD.values()]
    desc = dict(spearman_r_err_new=spearman([R[p] for p in S], [S[p]["err"] for p in S]) if len(S) > 2 else None,
                spearman_r_err_all=spearman(xs, ys))
    out = dict(profiles=S, verdicts=V, descriptive=desc)
    (ROOT / "analysis" / "bound.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"{'prof':4s} {'r':>5s} {'kU':>6s} {'K':>6s} {'err':>6s} {'A_cur':>6s} {'1-K':>6s} {'A+K':>5s} n unrel")
    for p in S:
        s = S[p]; print(f"{p:4s} {R[p]:5.2f} {KU[p]:6.3f} {s['kappa_median']:6.3f} {s['err']:+6.3f} {s['A_current']:6.2f} "
                        f"{s['predicted_A']:6.2f} {s['A_plus_K']:5.2f} {s['n']} {s['unreliable_share']:.2f}")
    print(json.dumps(V), json.dumps(desc))
