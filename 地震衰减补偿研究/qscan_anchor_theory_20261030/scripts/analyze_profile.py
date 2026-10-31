"""Score the Q-profile test (PLAN_PROFILE_ZH.md) from results7/ -> analysis/profile.json (+ printed summary)."""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
KAPPA = dict(S=0.880, D=0.294); CURV = dict(S=50.3, D=1.9); BASE_K = 0.675
RNG = np.random.default_rng(20261102); NB = 2000
VGRID = 1.0 / np.geomspace(0.4, 2.5, 15)   # 1/lib.QGRID


def load():
    recs = [json.loads(Path(f).read_text(encoding="utf-8")) for f in sorted(glob.glob(str(ROOT / "results7" / "*.json")))]
    return recs


def fit3(rows, key="v_sel"):
    X = np.array([[r["u"], r["v_anchor"], 1.0] for r in rows]); y = np.array([r[key] for r in rows])
    return np.linalg.lstsq(X, y, rcond=None)[0]


def boot(ids, fn):
    vals = []
    for _ in range(NB):
        pick = [ids[i] for i in RNG.integers(0, len(ids), len(ids))]; vals.append(fn(pick))
    vals = np.array(vals); return np.percentile(vals, 5, axis=0).tolist(), np.percentile(vals, 95, axis=0).tolist()


def summary(prof, recs):
    anch = [r for r in recs if r["profile"] == prof and r["scan"] != "w"]
    ids = sorted({r["id"] for r in anch}); by = {i: [r for r in anch if r["id"] == i] for i in ids}
    rows = lambda pick: [r for i in pick for r in by[i]]
    A, K, C = fit3(anch); lo, hi = boot(ids, lambda p: fit3(rows(p)))
    Ag, Kg, Cg = fit3(anch, "v_grid")
    one = [r for r in anch if r["scan"] == "a1.0"]
    s1 = float(np.polyfit([r["u"] for r in one], [r["v_sel"] for r in one], 1)[0])
    s1_ci = boot(ids, lambda p: float(np.polyfit([r["u"] for i in p for r in by[i] if r["scan"] == "a1.0"],
                                                [r["v_sel"] for i in p for r in by[i] if r["scan"] == "a1.0"], 1)[0]))
    bad = sum(r["edge"] or not r["neighbourhood_converged"] for r in anch) / len(anch)
    W = {r["id"]: r for r in recs if r["profile"] == prof and r["scan"] == "w"}

    def H(r):
        if any(x is None for x in r["J"]):
            return np.nan
        J = np.array(r["J"], float); jm = int(np.argmin(J)); lo_, hi_ = max(0, jm - 2), min(len(J), jm + 3)
        if hi_ - lo_ < 4:
            return np.nan
        return 2 * np.polyfit(VGRID[lo_:hi_], J[lo_:hi_], 2)[0]

    ratios = [H(r) / H(W[r["id"]]) for r in one if r["id"] in W]   # addendum: QScan-W dropped, descriptive only
    wsl = float(np.polyfit([r["u"] for r in W.values()], [r["v_sel"] for r in W.values()], 1)[0]) if W else None
    return dict(n_datasets=len(ids), n_scans=len(anch), kappa_pred=KAPPA[prof], A=float(A), K=float(K), C=float(C),
                A_plus_K=float(A + K), ci_A_K_C=[lo, hi], A_plus_K_ci=[lo[0] + lo[1], hi[0] + hi[1]],
                grid=dict(A=float(Ag), K=float(Kg), C=float(Cg)), slope_current=s1, slope_current_ci=list(s1_ci),
                unreliable_share=bad, curvature_ratio_median=float(np.nanmedian(ratios)) if ratios else None,
                curvature_ratios=[float(x) for x in ratios], curvature_pred=CURV[prof], qscanw_slope=wsl,
                per_dataset={i: {r["scan"]: r["v_sel"] for r in by[i]} | {"u": by[i][0]["u"]} for i in ids})


def verdicts(S):
    V = {}
    s, d = S.get("S"), S.get("D")
    ok = lambda x: x is not None and x["unreliable_share"] <= 0.25
    if ok(s) and ok(d):
        V["Q1"] = bool(abs(s["K"] - 0.880) <= 0.12 and abs(d["K"] - 0.294) <= 0.12)
        V["Q2"] = bool(d["K"] < BASE_K < s["K"] and d["ci_A_K_C"][1][1] < BASE_K and s["ci_A_K_C"][0][1] > BASE_K)
        V["Q3"] = bool(abs(s["A_plus_K"] - 1) <= 0.2 and abs(d["A_plus_K"] - 1) <= 0.2)
        V["Q4"] = bool(0.5 <= d["slope_current"] <= 0.9 and -0.05 <= s["slope_current"] <= 0.30)
    # Q5 cancelled by the addendum (QScan-W scans dropped for run time)
    return V


if __name__ == "__main__":
    R = load(); S = {p: summary(p, R) for p in ("S", "D") if any(r["profile"] == p for r in R)}
    S["verdicts"] = verdicts(S)
    (ROOT / "analysis" / "profile.json").write_text(json.dumps(S, indent=1), encoding="utf-8")
    for p in ("S", "D"):
        if p in S:
            print(p, {k: v for k, v in S[p].items() if k not in ("per_dataset", "curvature_ratios")})
    print("verdicts", S["verdicts"])
