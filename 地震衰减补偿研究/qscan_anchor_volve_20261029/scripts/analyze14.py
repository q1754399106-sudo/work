"""Experiment 14 analysis (PLAN_ZH.md section 4) -> analysis/calibration.json, analysis/summary.json, printed verdicts.

Calibration (development set only): offset(lam) = mean log(c_sel / c_true) of QScan-W over the 20 sections of
experiment 12; offset_SR from the SR factors of experiment 12 on the same sections (results/posthoc_qscan, git HEAD).
Calibrated estimates on the test sets: c_sel * exp(-offset).
Bootstrap: 2000 resamples of datasets within the group, seed 20261030, 90 % percentile intervals.
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
E12 = "地震衰减补偿研究/baselines_volve_20261029/results/posthoc_qscan"
NB, SEED, MARGIN, H1_MIN = 2000, 20261030, 0.02, 0.2


def load(stage):
    return [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / stage).glob("*.json"))]


def slope(sel, ct):
    return float(np.polyfit(np.log(ct), np.log(sel), 1)[0])


def boot(fn, n):
    rng = np.random.default_rng(SEED); vals = []
    for _ in range(NB):
        idx = rng.integers(0, n, n)
        if np.ptp(idx) == 0:
            continue
        vals.append(fn(idx))
    vals = np.array(vals); vals = vals[np.isfinite(vals)]
    return [float(np.percentile(vals, 5)), float(np.percentile(vals, 95))]


def verdict(ci, point, lo_rule, hi_rule):
    if lo_rule(ci, point):
        return "支持"
    if hi_rule(ci, point):
        return "相反"
    return "无法判定"


def describe(sel, ct):
    e = np.log(sel / ct)
    return dict(slope=slope(sel, ct), bias_log=float(e.mean()), abs_log=float(np.abs(e).mean()),
                within25=float(np.mean(np.abs(e) <= np.log(1.25))))


dev = load("dev"); ct_d = np.array([r["c_true"] for r in dev]); lams = dev[0]["lams"]
e12 = [json.loads(subprocess.run(["git", "show", f"HEAD:{E12}/{r['id']}.json"], cwd=REPO, capture_output=True,
                                 check=True).stdout.decode("utf-8")) for r in dev]
cal = {str(l): float(np.mean(np.log(np.array([r["qscan_w"][str(l)] for r in dev]) / ct_d))) for l in lams}
cal["SR"] = float(np.mean(np.log(np.array([d["c_sr"] for d in e12]) / ct_d)))
cal["dev_describe"] = {str(l): describe(np.array([r["qscan_w"][str(l)] for r in dev]), ct_d) for l in lams}
cal["dev_converged"] = {str(l): float(np.mean([np.mean(r["diag"]["converged"][i]) for r in dev])) for i, l in enumerate(lams)}
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "calibration.json").write_text(json.dumps(cal, indent=1), encoding="utf-8")

test = load("test"); out = dict(calibration=cal, groups={}, verdicts={})
for grp, members in (("T", ("T",)), ("S", ("L", "N"))):
    rs = [r for r in test if r["group"] in members]; ct = np.array([r["c_true"] for r in rs]); n = len(rs)
    est = {"QScan": np.array([r["qscan"] for r in rs]), "SR": np.array([r["c_sr"] for r in rs]),
           "SR_cal": np.array([r["c_sr"] for r in rs]) * np.exp(-cal["SR"])}
    for l in lams:
        est[f"W{l}"] = np.array([r["qscan_w"][str(l)] for r in rs])
        est[f"W{l}_cal"] = est[f"W{l}"] * np.exp(-cal[str(l)])
    g = dict(n=n, ids=[r["id"] for r in rs], describe={k: describe(v, ct) for k, v in est.items()},
             converged={str(l): float(np.mean([np.mean(r["diag"]["converged"][i]) for r in rs])) for i, l in enumerate(lams)},
             estimates={k: v.tolist() for k, v in est.items()})
    err = {k: np.abs(np.log(v / ct)) for k, v in est.items()}
    d3 = err[f"W{2 ** -4}_cal"] - err["SR_cal"]; d3r = err[f"W{2 ** -4}_cal"] - err["SR"]
    g["H3_diff_vs_SR_cal"] = dict(mean=float(d3.mean()), ci=boot(lambda i: d3[i].mean(), n))
    g["H3_diff_vs_SR_raw"] = dict(mean=float(d3r.mean()), ci=boot(lambda i: d3r[i].mean(), n))
    out["groups"][grp] = g
    ni = lambda ci, p: ci[1] <= MARGIN
    inf = lambda ci, p: ci[0] > MARGIN
    key = "H3" if grp == "T" else "H4"
    out["verdicts"][key] = dict(diff=g["H3_diff_vs_SR_cal"], result=verdict(g["H3_diff_vs_SR_cal"]["ci"], None, ni, inf))
    if grp == "T":
        a, b = est[f"W{0.5}"], est["QScan"]
        dS = slope(a, ct) - slope(b, ct); ci = boot(lambda i: slope(a[i], ct[i]) - slope(b[i], ct[i]), n)
        out["verdicts"]["H1"] = dict(diff=dS, ci=ci, result=verdict(ci, dS, lambda c, p: c[0] > 0 and p >= H1_MIN,
                                                                     lambda c, p: c[1] < 0))
        bias = [g["describe"][f"W{l}"]["bias_log"] for l in lams]; slp = [g["describe"][f"W{l}"]["slope"] for l in lams]
        hi, lo = est[f"W{lams[-1]}"], est[f"W{lams[0]}"]
        db = float(np.mean(np.log(hi / ct)) - np.mean(np.log(lo / ct)))
        cb = boot(lambda i: float(np.mean(np.log(hi[i] / ct[i])) - np.mean(np.log(lo[i] / ct[i]))), n)
        ds = slope(hi, ct) - slope(lo, ct); cs = boot(lambda i: slope(hi[i], ct[i]) - slope(lo[i], ct[i]), n)
        mono = all(np.diff(bias) > 0) and all(np.diff(slp) > 0)
        res = "支持" if mono and cb[0] > 0 and cs[0] > 0 else ("相反" if cb[1] < 0 or cs[1] < 0 else "无法判定")
        out["verdicts"]["H2"] = dict(bias=bias, slope=slp, monotone=mono, bias_diff=db, bias_ci=cb, slope_diff=ds,
                                     slope_ci=cs, result=res)
(ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print("calibration:", {k: round(v, 3) for k, v in cal.items() if not isinstance(v, dict)})
for grp, g in out["groups"].items():
    print(f"group {grp} (n={g['n']}), converged {g['converged']}")
    for k, s in g["describe"].items():
        print(f"  {k:<14} slope={s['slope']:5.2f} bias={s['bias_log']:+.3f} |err|={s['abs_log']:.3f} <=25%={s['within25']:.2f}")
for k, v in out["verdicts"].items():
    print(k, v["result"], {kk: vv for kk, vv in v.items() if kk != "result"})
