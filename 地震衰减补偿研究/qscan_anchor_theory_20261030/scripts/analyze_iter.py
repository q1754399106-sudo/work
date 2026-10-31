"""Score the iterated re-anchoring test (PLAN_ITER_ZH.md) from results5/ -> analysis/iter.json (+ printed summary)."""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
KAPPA = dict(main=0.675, mid=0.964, deep=1.140, main16=0.675)
RNG = np.random.default_rng(20261101); NB = 2000


def load():
    recs = {}
    for f in sorted(glob.glob(str(ROOT / "results5" / "*.json"))):
        r = json.loads(Path(f).read_text(encoding="utf-8")); recs.setdefault(r["arm"], []).append(r)
    return recs


def path(r):
    return np.array([1.0] + [h["v_sel"] for h in r["history"]])


def ratio(recs):
    num = den = 0.0
    for r in recs:
        d = np.diff(path(r)); num += float(np.sum(d[1:] * d[:-1])); den += float(np.sum(d[:-1] ** 2))
    return num / den if den > 0 else float("nan")


def boot(recs, fn):
    n = len(recs); vals = [fn([recs[i] for i in RNG.integers(0, n, n)]) for _ in range(NB)]
    return [float(np.nanpercentile(vals, 5)), float(np.nanpercentile(vals, 95))]


def fit(recs, which):
    u = np.array([r["u"] for r in recs]); v = np.array([which(r) for r in recs]); A, C = np.polyfit(u, v, 1)
    return float(A), float(C)


def monotone_conv(r):
    d = np.diff(path(r)); s = np.sign(d[0]); return bool(all(np.sign(x) == s for x in d if abs(x) > 0.01))


def monotone_away(r):
    dev = np.abs(path(r) - 1.0); return bool(np.all(np.diff(dev) > 0))


def unreliable(recs):
    steps = [h for r in recs for h in r["history"]]
    return sum(not h["neighbourhood_converged"] for h in steps) / max(1, len(steps))


def extrap(r, k):
    return (r["history"][0]["v_sel"] - k) / (1 - k)


def arm_summary(arm, recs):
    k = KAPPA[arm]; out = dict(n=len(recs), kappa=k, r=ratio(recs), r_ci=boot(recs, ratio),
                              unreliable_share=unreliable(recs), steps=[len(r["history"]) for r in recs],
                              edge_hits=sum(r["history"][-1]["edge"] for r in recs))
    out["monotone_converging"] = sum(monotone_conv(r) for r in recs)
    out["monotone_away"] = sum(monotone_away(r) for r in recs)
    if arm in ("main", "main16"):
        last = lambda r: r["history"][-1]["v_sel"]
        out["A_star"], out["C_line"] = fit(recs, last); out["A_star_ci"] = boot(recs, lambda s: fit(s, last)[0])
        out["A_step1"], out["C_step1"] = fit(recs, lambda r: r["history"][0]["v_sel"])
        cstar = lambda s: float(np.mean([last(r) - r["u"] for r in s]))
        out["C_star"] = cstar(recs); out["C_star_ci"] = boot(recs, cstar)
        dev = np.array([last(r) - extrap(r, k) for r in recs])
        out["extrap_absdev_median"] = float(np.median(np.abs(dev)))
        out["slope_vfinal_on_extrap"] = float(np.polyfit([extrap(r, k) for r in recs], [last(r) for r in recs], 1)[0])
        out["beta_step1_mean"] = float(np.mean([r["history"][0]["v_sel"] - (1 - k) * r["u"] - k for r in recs]))
        err = lambda v, u: abs(np.log(u / v))      # |log c error| (c = 1/v)
        out["median_abs_log_err"] = dict(
            step1=float(np.median([err(r["history"][0]["v_sel"], r["u"]) for r in recs])),
            final=float(np.median([err(last(r), r["u"]) for r in recs])),
            extrap=float(np.median([err(extrap(r, k), r["u"]) for r in recs])))
        out["by_c"] = {str(c): float(np.mean([last(r) - r["u"] for r in recs if r["c_true"] == c]))
                       for c in sorted({r["c_true"] for r in recs})}
        out["per_dataset"] = {r["id"]: dict(u=r["u"], path=path(r).tolist(), extrap=extrap(r, k)) for r in recs}
    else:
        out["per_dataset"] = {r["id"]: dict(u=r["u"], path=path(r).tolist(), edge=r["history"][-1]["edge"])
                              for r in recs}
    return out


def qscanw(recs):
    """Descriptive: QScan-W (lambda 1/16, T1 = 3.0) selections of the analysis-window test on seeds 9941-9944."""
    tr = json.loads((ROOT / "analysis" / "trunc.json").read_text(encoding="utf-8"))["qscanw"]["3.0"]["sel"]
    w = {f"c{c}_s{s}": 1 / cs for c, s, cs in tr}; rows = {}
    for r in recs:
        if r["id"] in w:
            rows[r["id"]] = dict(u=r["u"], qscanw=w[r["id"]], final=r["history"][-1]["v_sel"],
                                 step1=r["history"][0]["v_sel"], extrap=extrap(r, KAPPA["main"]))
    med = {k: float(np.median([abs(np.log(v["u"] / v[k])) for v in rows.values()]))
           for k in ("qscanw", "final", "step1", "extrap")} if rows else {}
    return dict(rows=rows, median_abs_log_err=med)


def verdicts(S):
    m, md, d = S.get("main"), S.get("mid"), S.get("deep"); V = {}
    ok = lambda s: s is not None and s["unreliable_share"] <= 0.25
    if ok(m):
        V["P1"] = bool(0.55 <= m["r"] <= 0.80 and m["monotone_converging"] >= 20)
        V["P3"] = bool(0.8 <= m["A_star"] <= 1.2)
        V["P4"] = bool(m["extrap_absdev_median"] <= 0.08 and 0.8 <= m["slope_vfinal_on_extrap"] <= 1.2)
        V["P5"] = bool(m["C_star_ci"][0] > 0 and 0.05 <= m["C_star"] <= 0.5)
    if ok(m) and ok(md) and d is not None:
        V["P2"] = bool(0.84 <= md["r"] <= 1.08 and 1.02 <= d["r"] <= 1.26 and m["r"] < md["r"] < d["r"]
                       and d["monotone_away"] >= round(0.75 * d["n"]))   # addendum: 5/6
    return V


if __name__ == "__main__":
    R = load(); S = {arm: arm_summary(arm, recs) for arm, recs in R.items()}
    if "main" in R:
        S["qscanw_compare"] = qscanw(R["main"])
    S["verdicts"] = verdicts(S)
    (ROOT / "analysis" / "iter.json").write_text(json.dumps(S, indent=1), encoding="utf-8")
    for arm, s in S.items():
        if arm in KAPPA:
            print(arm, {k: v for k, v in s.items() if k not in ("per_dataset",)})
    print("qscanw", S.get("qscanw_compare", {}).get("median_abs_log_err")); print("verdicts", S["verdicts"])
