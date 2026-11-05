"""Score PLAN_ZH.md from ../results -> ../analysis/orig.json."""
import json, glob, sys
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]; RNG = np.random.default_rng(20261106); NB = 2000
KU = {"W1": 0.075, "W2": 0.375, "W3": 0.75}; DVA = 0.01


def num(x):
    return float(x) if isinstance(x, (int, float)) else np.nan


def load():
    T = {}
    for f in glob.glob(str(ROOT / "results" / "*.json")):
        r = json.loads(Path(f).read_text(encoding="utf-8"))
        key = "known" if r["cond"] == "known" else f"{r['cond']}_{num(r['va']):.2f}"
        T[(key, r["Qt"], r["seed"])] = r
    return T


def slope(T, key, crit, qts, seeds):
    x, y = [], []
    for s in seeds:
        xs = np.array([1 / q for q in qts]); ys = np.array([num(T[(key, q, s)][crit]) for q in qts]); ok = np.isfinite(ys)
        if ok.sum() >= 2:
            x += list(xs[ok] - xs[ok].mean()); y += list(ys[ok] - ys[ok].mean())
    x, y = np.array(x), np.array(y)
    return float(x @ y / (x @ x)) if len(x) else np.nan


def kappa(T, w, crit, qts, seeds):
    d = np.array([(num(T[(f"{w}_0.01", q, s)][crit]) - num(T[(f"{w}_0.00", q, s)][crit])) / DVA for s in seeds for q in qts])
    return float(np.median(d[np.isfinite(d)])) if np.isfinite(d).any() else np.nan


def boot(fn, seeds):
    b = np.array([fn(list(RNG.choice(seeds, len(seeds)))) for _ in range(NB)])
    return [float(np.nanpercentile(b, 5)), float(np.nanpercentile(b, 95))]


if __name__ == "__main__":
    T = load(); qts = sorted({k[1] for k in T}); seeds = sorted({k[2] for k in T}); keys = sorted({k[0] for k in T})
    out = dict(n=len(T), conds={})
    for crit in ("v1", "v2"):
        for key in keys:
            vals = np.array([num(T[(key, q, s)][crit]) for q in qts for s in seeds])
            e = dict(missing=float(np.mean(~np.isfinite(vals))), A=slope(T, key, crit, qts, seeds))
            e["A_ci"] = boot(lambda ss: slope(T, key, crit, qts, ss), seeds)
            e["median_Q_by_Qt"] = {str(q): float(np.nanmedian([1 / num(T[(key, q, s)][crit]) for s in seeds])) for q in qts}
            out["conds"][f"{crit}:{key}"] = e
        for w in KU:
            out["conds"][f"{crit}:K_{w}"] = dict(K=kappa(T, w, crit, qts, seeds), KU=KU[w],
                                                 K_ci=boot(lambda ss: kappa(T, w, crit, qts, ss), seeds))
    c = out["conds"]; ok = lambda *ks: all(c[f"v2:{k}"]["missing"] <= 0.25 for k in ks)
    dA = lambda ss: slope(T, "W1_0.00", "v2", qts, ss) - slope(T, "W3_0.00", "v2", qts, ss)
    d = dA(seeds); dci = boot(dA, seeds); V = {}
    V["O1"] = bool(0.8 <= c["v2:known"]["A"] <= 1.2) if ok("known") else None
    V["O2"] = dict(diff=d, ci=dci, holds=bool(d >= 0.2 and dci[0] > 0)) if ok("W1_0.00", "W3_0.00") else None
    V["O3"] = dict(W1=bool(abs(c["v2:K_W1"]["K"] - 0.075) <= 0.15), W2=bool(abs(c["v2:K_W2"]["K"] - 0.375) <= 0.2),
                   W3=bool(abs(c["v2:K_W3"]["K"] - 0.75) <= 0.2))
    V["O4"] = bool(abs(c["v2:W3_0.00"]["A"] - 0.25) <= 0.2) if ok("W3_0.00") else None
    k = c["v1:K_W3"]; lo, hi = k["K_ci"]
    V["C1_class"] = ("anchored" if (k["K"] >= 0.3 and lo > 0) else
                     "cancels" if (abs(k["K"]) <= 0.15 and lo >= -0.3 and hi <= 0.3) else "undecided")
    out["verdicts"] = V; (ROOT / "analysis" / "orig.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for crit, nm in (("v1", "C1"), ("v2", "C2")):
        print(nm)
        for key in keys:
            e = c[f"{crit}:{key}"]
            print(f"  {key:9s} A {e['A']:+.3f} [{e['A_ci'][0]:+.2f},{e['A_ci'][1]:+.2f}] missing {e['missing']:.2f} "
                  f"medQ {[round(v) for v in e['median_Q_by_Qt'].values()]}")
        for w in KU:
            e = c[f"{crit}:K_{w}"]; print(f"  K_{w} {e['K']:+.3f} [{e['K_ci'][0]:+.2f},{e['K_ci'][1]:+.2f}] (kU {e['KU']})")
    print(json.dumps(V))
