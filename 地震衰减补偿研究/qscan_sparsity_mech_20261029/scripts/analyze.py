"""Summary of the three prototypes -> analysis/summary.json, printed table.

proto/            reflectivity family (sparse / log / gauss) x QScan weights 0.5, 2^-5, plus SR (6 datasets each)
consistent/       QScan with the wavelet re-estimated per factor, weights 0.5, 0.125 (20 sections of experiment 12)
consistent_lams/  the same at weights 2^-4, 2^-6
References for the 20 sections: QScan and SR of experiment 12 (results/posthoc_qscan, git HEAD).
Per row: log-log slope on the true factor, mean signed and absolute log error, share within 25 %, leave-one-out
offset-corrected |log error|, converged share, boundary share.
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
E12 = "地震衰减补偿研究/baselines_volve_20261029/results/posthoc_qscan"


def load(d):
    return [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / d).glob("*.json"))]


def stats(sel, ct, conv=None):
    e = np.log(sel / ct); loo = np.array([e[k] - np.delete(e, k).mean() for k in range(len(e))])
    return dict(n=len(e), slope=float(np.polyfit(np.log(ct), np.log(sel), 1)[0]), bias_log=float(e.mean()),
                abs_log=float(np.abs(e).mean()), within25=float(np.mean(np.abs(e) <= np.log(1.25))),
                loo_offset_abs_log=float(np.abs(loo).mean()), converged=conv,
                boundary=float(np.mean((sel <= 0.401) | (sel >= 2.49))))


rows = {}
P = load("proto")
for fam in ("sparse", "log", "gauss"):
    rs = [r for r in P if r["family"] == fam]; ct = np.array([r["c_true"] for r in rs])
    for i, lam in enumerate(rs[0]["lams"]):
        rows[f"proto {fam} QScan lam={lam}"] = stats(np.array([r["selected"][str(lam)] for r in rs]), ct,
                                                     float(np.mean([np.mean(r["converged"][i]) for r in rs])))
    rows[f"proto {fam} SR"] = stats(np.array([r["c_sr"] for r in rs]), ct)
for d in ("consistent", "consistent_lams"):
    rs = load(d); ct = np.array([r["c_true"] for r in rs])
    for i, lam in enumerate(rs[0]["lams"]):
        rows[f"consistent lam={lam}"] = stats(np.array([r["selected"][str(lam)] for r in rs]), ct,
                                              float(np.mean([np.mean(r["converged"][i]) for r in rs])))
e12 = [json.loads(subprocess.run(["git", "show", f"HEAD:{E12}/{r['id']}.json"], cwd=REPO, capture_output=True,
                                 check=True).stdout.decode("utf-8")) for r in rs]
rows["exp12 QScan lam=0.5"] = stats(np.array([d["qscan"] for d in e12]), ct)
rows["exp12 SR"] = stats(np.array([d["c_sr"] for d in e12]), ct)
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "summary.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
print(f"{'row':<34} {'n':>3} {'slope':>5} {'bias':>7} {'|err|':>6} {'<=25%':>5} {'LOO':>6} {'conv':>5}")
for k, s in rows.items():
    cv = "  -  " if s["converged"] is None else f"{s['converged']:5.2f}"
    print(f"{k:<34} {s['n']:3d} {s['slope']:5.2f} {s['bias_log']:+7.3f} {s['abs_log']:6.3f} {s['within25']:5.2f} "
          f"{s['loo_offset_abs_log']:6.3f} {cv}")
