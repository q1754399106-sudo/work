"""Verdicts of PLAN_CONFIRM_S_ZH.md from results/confirmS -> analysis/confirmS.json and printed tables.
Per state, scores are first averaged over the 9 Q multipliers (equal weight) or over a side (under: m 1.4, 1.7, 2.0).
Semi: unit = dataset (12), 90% interval from 2000 bootstrap resamples of datasets (seed 20261199). FD: unit = model (6), the 3
realisations averaged first; mean and the number of models with a positive difference.
usage: analyze_confirmS.py [--dev DIR]   (--dev: read another results directory, for the code check on development output)"""
import sys, json, glob
import numpy as np
from env import ROOT

MS = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); UNDER = (1.4, 1.7, 2.0); OVER = (0.5, 0.6, 0.7); NEAR = (0.85, 1.0, 1.2)
SRC = sys.argv[sys.argv.index("--dev") + 1] if "--dev" in sys.argv else str(ROOT / "results" / "confirmS")


def table(prefix, key="S"):
    """{(var, state, method): array (units, 9)}; FD: realisations averaged per model."""
    out = {}
    if prefix == "S6":
        files = sorted(glob.glob(f"{SRC}/S6_*.json"))
        for i, p in enumerate(files):
            for r in json.load(open(p))["records"]:
                out.setdefault((r["var"], r["state"], r["method"]), np.full((len(files), 9), np.nan))[i, MS.index(r["m"])] = r[key]
        return out, [p.split("S6_")[1][:-5] for p in files]
    files = sorted(glob.glob(f"{SRC}/FD3_*.json")); models = sorted({p.split("FD3_")[1].split("_")[0] for p in files}, key=lambda h: int(h[1:]))
    acc = {}
    for p in files:
        h, lv = p.split("FD3_")[1].split("_")[:2]
        for r in json.load(open(p))["records"]:
            acc.setdefault((lv, r["method"]), {}).setdefault(h, np.zeros(9)); acc[(lv, r["method"])][h][MS.index(r["m"])] += r[key] / 3
    return {("FD", lv, meth): np.array([v[h] for h in models]) for (lv, meth), v in acc.items()}, models


def side(a, ms):
    return a[:, [MS.index(m) for m in ms]].mean(1)


def boot(d):
    g = np.random.default_rng(20261199); return [float(x) for x in np.percentile(d[g.integers(0, len(d), (2000, len(d)))].mean(1), [5, 95])]


def diff(T, var, st, a, b, ms=MS):
    d = side(T[(var, st, a)], ms) - side(T[(var, st, b)], ms); return dict(mean=float(d.mean()), ci=boot(d), pos=int((d > 0).sum()), n=len(d))


S6, ids = table("S6"); FD, models = table("FD3"); T = {**S6, **FD}; res = dict(datasets=ids, models=models, rules={}, describe={})
R = res["rules"]
R["U1"] = dict(test="semi cal: UQWS - UQW, 9 m", **diff(T, "O", "cal", "UQWS", "UQW")); R["U1"]["holds"] = R["U1"]["ci"][0] > 0
for st in ("res4", "gau4"):
    R["U2_" + st] = dict(test=f"semi {st}: UQWS - UQW, 9 m", **diff(T, "O", st, "UQWS", "UQW")); R["U2_" + st]["holds"] = R["U2_" + st]["ci"][0] > -0.01
for lv in ("30dB", "14dB"):
    R["U3_" + lv] = dict(test=f"FD {lv}: UQWS - UQW, 9 m", **diff(T, "FD", lv, "UQWS", "UQW")); R["U3_" + lv]["holds"] = R["U3_" + lv]["mean"] >= -0.01
R["U4"] = dict(test="semi cal, m 1.4-2.0: UQWS - UQW", **diff(T, "O", "cal", "UQWS", "UQW", UNDER)); R["U4"]["holds"] = R["U4"]["ci"][0] > 0
R["U5_semi"] = dict(test="semi res4: UQWS - VDBin, 9 m", **diff(T, "O", "res4", "UQWS", "VDBin")); R["U5_semi"]["holds"] = R["U5_semi"]["ci"][0] > 0
R["U5_fd"] = dict(test="FD 14dB: UQWS - VDBin, 9 m", **diff(T, "FD", "14dB", "UQWS", "VDBin")); R["U5_fd"]["holds"] = R["U5_fd"]["mean"] > 0 and R["U5_fd"]["pos"] >= 5
res["confirmed"] = all(v["holds"] for v in R.values())
states = [("O", s) for s in ("cal", "res4", "gau4")] + [("F", s) for s in ("cal", "res4", "gau4")] + [("FD", "30dB"), ("FD", "14dB")]
for var, st in states:
    if (var, st, "UQW") not in T:
        continue
    D = res["describe"][f"{var}_{st}"] = {}
    D["mean_by_m"] = {k: [float(x) for x in T[(var, st, k)].mean(0)] for (v_, s_, k) in T if v_ == var and s_ == st}
    for a, b in (("UQWS", "UQW"), ("UQWS", "VDBin"), ("UQWS", "VDB"), ("UQW", "VDBin"), ("UQW", "VDB"), ("VDBin", "VDB"), ("UQWS", "HWin"), ("UQWS", "IQ70")):
        D[f"{a}-{b}"] = {nm: diff(T, var, st, a, b, ms) for nm, ms in (("all", MS), ("over", OVER), ("near", NEAR), ("under", UNDER))}
    if var != "FD":
        D["T2_UQWS-UQW"] = (lambda d: dict(mean=float(d.mean()), pos=int((d > 0).sum())))(side(table("S6", "T2")[0][(var, st, "UQWS")], MS) - side(table("S6", "T2")[0][(var, st, "UQW")], MS))
if glob.glob(f"{SRC}/S6_*.json"):
    res["bands"] = [dict(id=json.load(open(p))["id"], bands=json.load(open(p))["bands"]) for p in sorted(glob.glob(f"{SRC}/S6_*.json"))]
if "--dev" not in sys.argv:
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "confirmS.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
for k, v in R.items():
    print("%-9s %-36s %+.3f  [%+.3f, %+.3f]  %d/%d  %s" % (k, v["test"], v["mean"], *v["ci"], v["pos"], v["n"], "holds" if v["holds"] else "FAILS"))
print("confirmed:", res["confirmed"])
for var, st in states:
    D = res["describe"].get(f"{var}_{st}")
    if D:
        print("%s %-5s " % (var, st) + "  ".join("%s %+.3f" % (k, D[k]["all"]["mean"]) for k in ("UQWS-UQW", "UQWS-VDBin", "UQW-VDBin", "UQW-VDB", "VDBin-VDB", "UQWS-HWin")))
