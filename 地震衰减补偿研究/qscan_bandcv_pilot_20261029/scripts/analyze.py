"""Summary of the band cross-validation pilot -> analysis/summary.json, printed table.

Per criterion (E_raw, E_gain) and weight: mean signed log error, mean |log error|, share within 25 %, Spearman rank
correlation, log-log slope of selected on true c, boundary share, converged share; plus leave-one-out offset
correction (subtract the mean log error of the other 19) since a constant bias is calibratable. Reference rows:
QScan (full-band lasso, weight 0.5) and SR, both read from experiment 12 (results/posthoc_qscan, git HEAD).
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
E12 = "地震衰减补偿研究/baselines_volve_20261029/results/posthoc_qscan"
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results").glob("T_*.json"))]
lams = recs[0]["lams"]; grid = np.array(recs[0]["grid"]); ct = np.array([r["c_true"] for r in recs])
e12 = [json.loads(subprocess.run(["git", "show", f"HEAD:{E12}/{r['id']}.json"], cwd=REPO, capture_output=True,
                                 check=True).stdout.decode("utf-8")) for r in recs]


def rank(x):
    return np.argsort(np.argsort(x, kind="stable"), kind="stable").astype(float)


def stats(sel, boundary=True):
    e = np.log(sel / ct); x = np.log(ct)
    loo = np.array([e[k] - np.delete(e, k).mean() for k in range(len(e))])
    return dict(bias_log=float(e.mean()), abs_log=float(np.abs(e).mean()), within25=float(np.mean(np.abs(e) <= np.log(1.25))),
                spearman=float(np.corrcoef(rank(sel), rank(ct))[0, 1]), slope=float(np.polyfit(x, np.log(sel), 1)[0]),
                resid_sd=float(np.std(e - e.mean(), ddof=1)),
                loo_offset_abs_log=float(np.abs(loo).mean()), loo_offset_within25=float(np.mean(np.abs(loo) <= np.log(1.25))),
                boundary=float(np.mean((sel == grid[0]) | (sel == grid[-1]))) if boundary else None,
                selected=dict(zip([r["id"] for r in recs], sel.tolist())))


out = dict(n=len(recs), rows={})
for m in ("E_raw", "E_gain"):
    for i, lam in enumerate(lams):
        s = stats(np.array([r["selected"][m][str(lam)] for r in recs]))
        s["converged"] = float(np.mean([np.mean(r["converged"][i]) for r in recs])); out["rows"][f"{m} lam={lam}"] = s
out["rows"]["QScan (exp 12)"] = stats(np.array([d["qscan"] for d in e12]))
out["rows"]["SR (exp 12)"] = stats(np.array([d["c_sr"] for d in e12]), boundary=False)
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"{'row':<22} {'bias':>7} {'|err|':>6} {'<=25%':>6} {'rho':>5} {'slope':>5} {'sd':>5} {'LOO|err|':>8} {'LOO25':>5} {'bound':>5}")
for k, s in out["rows"].items():
    b = "  -  " if s["boundary"] is None else f"{s['boundary']:5.2f}"
    print(f"{k:<22} {s['bias_log']:+7.3f} {s['abs_log']:6.3f} {s['within25']:6.2f} {s['spearman']:5.2f} {s['slope']:5.2f} "
          f"{s['resid_sd']:5.3f} {s['loo_offset_abs_log']:8.3f} {s['loo_offset_within25']:5.2f} {b}")
