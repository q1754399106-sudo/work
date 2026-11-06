"""Score PLAN_ZH.md from ../results -> ../analysis/factor.json."""
import json, glob, sys
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]; RNG = np.random.default_rng(20261107); NB = 2000; DVA = 0.01; F0 = 30.0
CONDS = ("known", "W3_0", "W3_1", "F0.8", "F1.25"); NM = dict(vJ="J", v1="E17", vd="chengD", v2="C2")


def num(x):
    return float(x) if isinstance(x, (int, float)) else np.nan


def load():
    T = {}
    for f in glob.glob(str(ROOT / "results" / "*.json")):
        r = json.loads(Path(f).read_text(encoding="utf-8"))
        if "solver" in r:
            T[(r["solver"], r["Qt"], r["seed"], r["cond"])] = r
    return T


def K(T, sv, cr, qts, seeds):
    d = np.array([(num(T[(sv, q, s, "W3_1")][cr]) - num(T[(sv, q, s, "W3_0")][cr])) / DVA for s in seeds for q in qts])
    return float(np.median(d[np.isfinite(d)])) if np.isfinite(d).any() else np.nan


def slope(T, sv, cr, cond, qts, seeds):
    x, y = [], []
    for s in seeds:
        xs = np.array([1 / q for q in qts]); ys = np.array([num(T[(sv, q, s, cond)][cr]) for q in qts]); ok = np.isfinite(ys)
        if ok.sum() >= 2:
            x += list(xs[ok] - xs[ok].mean()); y += list(ys[ok] - ys[ok].mean())
    x, y = np.array(x), np.array(y)
    return float(x @ y / (x @ x)) if len(x) else np.nan


def ratio_slope(T, sv, cr, qts, seeds, cross_only):
    """log(Q_sel/Q_true) on log(fw/f0), fixed effect per (Q_true, seed); a repeated seed counts again."""
    X, Y = [], []
    for s in seeds:
        for q in qts:
            x, y = [], []
            for c in CONDS:
                r = T[(sv, q, s, c)]; v = num(r[cr])
                if np.isfinite(v) and v > 0 and (r["cross"] or not cross_only):
                    x.append(np.log(r["fw"] / F0)); y.append(np.log((1 / v) / q))
            if len(x) >= 2:
                X += list(np.array(x) - np.mean(x)); Y += list(np.array(y) - np.mean(y))
    X, Y = np.array(X), np.array(Y)
    return float(X @ Y / (X @ X)) if len(X) >= 4 else np.nan


def boot(fn, seeds):
    b = np.array([fn(list(RNG.choice(seeds, len(seeds)))) for _ in range(NB)])
    return [float(np.nanpercentile(b, 5)), float(np.nanpercentile(b, 95))]


def repo():
    rows = {}
    for f in glob.glob(str(ROOT / "results" / "repo_*.json")):
        r = json.loads(Path(f).read_text(encoding="utf-8")); rows.setdefault(r["dataset"], {})[r["c_anchor"]] = r
    dv = 1 / 0.7 - 1 / 1.4; out = {}
    for d, a in sorted(rows.items()):
        out[d] = dict(chengD=(a[0.7]["v_chengD"] - a[1.4]["v_chengD"]) / dv, e17=(a[0.7]["v_e17"] - a[1.4]["v_e17"]) / dv,
                      cross=[a[0.7]["cross_e17"], a[1.4]["cross_e17"]], f_in=[a[0.7]["f_in"], a[1.4]["f_in"]])
    return out


if __name__ == "__main__":
    T = load(); qts = [50, 70, 100, 140]; seeds = list(range(1201, 1207)); out = dict(n=len(T)); V = {}
    CR = dict(lasso=("vJ", "v1", "vd", "v2"), irls=("v1", "vd", "v2"))
    for sv, crs in CR.items():
        e = {}
        for cr in crs:
            f = lambda ss: K(T, sv, cr, qts, ss)
            e[cr] = dict(K=f(seeds), K_ci=boot(f, seeds), A_W3_0=slope(T, sv, cr, "W3_0", qts, seeds), A_known=slope(T, sv, cr, "known", qts, seeds),
                         rs=ratio_slope(T, sv, cr, qts, seeds, False), rs_ci=boot(lambda ss: ratio_slope(T, sv, cr, qts, ss, False), seeds),
                         medQ={c: [float(np.nanmedian([1 / num(T[(sv, q, s, c)][cr]) for s in seeds])) for q in qts] for c in CONDS})
        c = e["v1"]; c["rs_cross"] = ratio_slope(T, sv, "v1", qts, seeds, True); c["rs_cross_ci"] = boot(lambda ss: ratio_slope(T, sv, "v1", qts, ss, True), seeds)
        e["cross"] = {x: {str(q): float(np.mean([T[(sv, q, s, x)]["cross"] for s in seeds])) for q in qts} for x in CONDS}
        e["cross_all"] = float(np.mean([T[(sv, q, s, x)]["cross"] for q in qts for s in seeds for x in CONDS]))
        fd = lambda ss: K(T, sv, "vd", qts, ss) - K(T, sv, "v1", qts, ss); dci = boot(fd, seeds); dk = fd(seeds)
        usable = bool(0.6 <= e["vd"]["A_known"] <= 1.4)
        V[f"D_{sv}"] = dict(dK=dk, ci=dci, K_e17=c["K"], K_chengD=e["vd"]["K"], A_known_chengD=e["vd"]["A_known"], chengD_usable=usable,
                            holds=bool(dk >= 0.3 and dci[0] > 0 and c["K"] <= 0.15) if usable else None)
        V[f"F_{sv}"] = dict(s=c["rs_cross"], ci=c["rs_cross_ci"], cross_all=e["cross_all"],
                            holds=None if e["cross_all"] < 0.5 else bool(0.7 <= c["rs_cross"] <= 1.3 and c["rs_cross_ci"][0] >= 0.5 and c["rs_cross_ci"][1] <= 1.5))
        out[sv] = e
    fs = lambda ss: K(T, "lasso", "v1", qts, ss) - K(T, "irls", "v1", qts, ss); sci = boot(fs, seeds)
    V["S_e17"] = dict(diff=fs(seeds), ci=sci, holds=bool(abs(fs(seeds)) <= 0.2 and sci[0] >= -0.4 and sci[1] <= 0.4))
    k = out["lasso"]["vJ"]; V["J_lasso"] = dict(K=k["K"], ci=k["K_ci"], A=k["A_W3_0"], A_plus_K=k["K"] + k["A_W3_0"], holds=bool(abs(k["K"] - 0.75) <= 0.2))
    rp = repo(); d = np.array([x["chengD"] - x["e17"] for x in rp.values()]); out["repo"] = rp
    V["D_repo"] = dict(median=float(np.median(d)), positive=int((d > 0).sum()), n=len(d), holds=bool(np.median(d) >= 0.3 and (d > 0).sum() >= 5))
    out["verdicts"] = V; (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "factor.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for sv in CR:
        print(sv, "cross_all %.2f" % out[sv]["cross_all"], json.dumps(out[sv]["cross"]["W3_0"]))
        for cr in CR[sv]:
            e = out[sv][cr]
            print(f"  {NM[cr]:7s} K {e['K']:+.3f} [{e['K_ci'][0]:+.2f},{e['K_ci'][1]:+.2f}] A(W3_0) {e['A_W3_0']:+.2f} A(known) {e['A_known']:+.2f} s {e['rs']:+.2f} [{e['rs_ci'][0]:+.2f},{e['rs_ci'][1]:+.2f}]")
            for c in CONDS:
                print(f"      {c:6s} medQ {[round(x) for x in e['medQ'][c]]}")
    for k, x in rp.items():
        print(k, "chengD %+.3f e17 %+.3f" % (x["chengD"], x["e17"]), x["cross"], [round(f, 1) for f in x["f_in"]])
    print(json.dumps(V, indent=0))
