"""Score the second profile test (PLAN_PROFILE2_ZH.md) from results10/ -> analysis/profile2.json."""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RNG = np.random.default_rng(20261105); NB = 2000
KU = {"D2": 0.489, "M": 0.595}; DV = 1 / 0.7 - 1 / 1.4


def summary(prof, recs):
    by = {}
    for r in recs:
        if r["profile"] == prof:
            by.setdefault(r["id"], {})[r["c_anchor"]] = r
    ids = sorted(i for i in by if len(by[i]) == 3)
    kap = {i: (by[i][0.7]["v_sel"] - by[i][1.4]["v_sel"]) / DV for i in ids}
    u = {i: by[i][1.0]["u"] for i in ids}

    def stats(pick):
        k = np.median([kap[i] for i in pick]); A = np.polyfit([u[i] for i in pick], [by[i][1.0]["v_sel"] for i in pick], 1)[0]
        X = np.array([[by[i][a]["u"], 1 / a, 1.0] for i in pick for a in (0.7, 1.0, 1.4)])
        y = np.array([by[i][a]["v_sel"] for i in pick for a in (0.7, 1.0, 1.4)])
        Af, Kf, _ = np.linalg.lstsq(X, y, rcond=None)[0]
        return float(k), float(A), float(Af + Kf)

    k, A, AK = stats(ids); B = np.array([stats([ids[j] for j in RNG.integers(0, len(ids), len(ids))]) for _ in range(NB)])
    ci = lambda c: [float(np.percentile(B[:, c], 5)), float(np.percentile(B[:, c], 95))]
    bad = sum(r["edge"] or not r["neighbourhood_converged"] for i in ids for r in by[i].values()) / max(1, 3 * len(ids))
    return dict(n=len(ids), kappa_uniform=KU[prof], kappa_median=k, kappa_ci=ci(0), A_current=A, A_current_ci=ci(1),
                predicted_A=1 - k, A_plus_K=AK, A_plus_K_ci=ci(2), unreliable_share=bad,
                per_dataset={i: dict(u=u[i], kappa=kap[i], v={str(a): by[i][a]["v_sel"] for a in (0.7, 1.0, 1.4)})
                             for i in ids})


if __name__ == "__main__":
    recs = [json.loads(Path(f).read_text()) for f in glob.glob(str(ROOT / "results10" / "*.json"))]
    S = {p: summary(p, recs) for p in KU if any(r["profile"] == p for r in recs)}; V = {}
    ok = lambda p: p in S and S[p]["unreliable_share"] <= 0.25
    if ok("D2") and ok("M"):
        V["G1"] = bool(all(abs(S[p]["A_current"] - S[p]["predicted_A"]) <= 0.15 for p in KU))
        V["G2"] = bool(S["D2"]["kappa_median"] - KU["D2"] >= 0.08 and abs(S["M"]["kappa_median"] - KU["M"]) <= 0.10)
        V["G3"] = bool(all(abs(S[p]["A_plus_K"] - 1) <= 0.2 for p in KU))
    S["verdicts"] = V
    (ROOT / "analysis" / "profile2.json").write_text(json.dumps(S, indent=1), encoding="utf-8")
    for p in KU:
        if p in S:
            print(p, {k: v for k, v in S[p].items() if k != "per_dataset"})
    print("verdicts", V)
