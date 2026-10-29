"""Summary of the pilot: QScan selection against the lasso weight -> analysis/summary.json, printed table.

Per weight: mean signed log error log(c_sel / c_true), mean |log error|, rank correlation with c_true, share at the
grid boundary, mean selection per true c, share of non-converged solves. Check: weight 0.5 must reproduce the
QScan factors of experiment 12 (results/posthoc_qscan, read from git HEAD).
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
E12 = "地震衰减补偿研究/baselines_volve_20261029/results/posthoc_qscan"
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results").glob("T_*.json"))]
lams = recs[0]["lams"]; grid = np.array(recs[0]["grid"]); ct = np.array([r["c_true"] for r in recs])


def rank(x):
    return np.argsort(np.argsort(x, kind="stable"), kind="stable").astype(float)


def spearman(a, b):
    ra, rb = rank(a), rank(b)
    return float(np.corrcoef(ra, rb)[0, 1]) if ra.std() > 0 and rb.std() > 0 else float("nan")


out = dict(n=len(recs), c_true_levels=sorted(set(ct.tolist())), by_lambda={})
for i, lam in enumerate(lams):
    sel = np.array([r["selected"][str(lam)] for r in recs]); e = np.log(sel / ct)
    conv = np.mean([np.mean(r["converged"][i]) for r in recs])
    out["by_lambda"][str(lam)] = dict(
        bias_log=float(e.mean()), abs_log=float(np.abs(e).mean()), within25=float(np.mean(np.abs(e) <= np.log(1.25))),
        spearman=spearman(sel, ct), boundary=float(np.mean((sel == grid[0]) | (sel == grid[-1]))), converged=float(conv),
        mean_sel_by_c={str(c): float(np.exp(np.log(sel[ct == c]).mean())) for c in sorted(set(ct.tolist()))},
        selected={r["id"]: s for r, s in zip(recs, sel.tolist())})
chk = {}
for r in recs:
    try:
        e12 = json.loads(subprocess.run(["git", "show", f"HEAD:{E12}/{r['id']}.json"], cwd=REPO, capture_output=True,
                                        check=True).stdout.decode("utf-8"))
        chk[r["id"]] = abs(e12["qscan"] - r["selected"]["0.5"])
    except subprocess.CalledProcessError:
        chk[r["id"]] = None
out["check_lambda05_vs_e12"] = dict(max_abs=max(v for v in chk.values() if v is not None), missing=sum(v is None for v in chk.values()))
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print("check vs exp 12 at 0.5:", out["check_lambda05_vs_e12"])
levels = out["c_true_levels"]
print(f"{'lambda':>9} {'bias':>7} {'|err|':>6} {'<=25%':>6} {'rho':>5} {'bound':>5} {'conv':>5}  " + " ".join(f"c={c:<5}" for c in levels))
for lam, s in out["by_lambda"].items():
    print(f"{float(lam):9.4f} {s['bias_log']:+7.3f} {s['abs_log']:6.3f} {s['within25']:6.2f} {s['spearman']:5.2f} "
          f"{s['boundary']:5.2f} {s['converged']:5.2f}  " + " ".join(f"{s['mean_sel_by_c'][str(c)]:<7.3f}" for c in levels))
