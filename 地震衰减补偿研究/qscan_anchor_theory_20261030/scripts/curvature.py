"""Curvature check (PLAN_CURV_ZH.md): curvature of the QScan objective in v = 1/c, fixed anchor vs re-estimated wavelet.

H_a: anchor experiment, c_a = 1 (results/<id>__a1.0.json); H_w: QScan-W of experiment 14, development set (git HEAD).
Quadratic fit in v on the minimum and two grid points on each side (skipped when the minimum is within two points of
the grid edge); curvature = 2 x the quadratic coefficient. -> analysis/curvature.json
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
QGRID = np.exp(np.linspace(np.log(0.4), np.log(2.5), 15)); V = 1 / QGRID
PRED = 8.0


def git_json(path):
    return json.loads(subprocess.run(["git", "show", f"HEAD:地震衰减补偿研究/{path}"], cwd=REPO, capture_output=True, check=True).stdout)


def curv(J):
    J = np.asarray(J, float); j = int(np.argmin(J))
    if j < 2 or j > len(J) - 3:
        return None
    return float(2 * np.polyfit(V[j - 2:j + 3], J[j - 2:j + 3], 2)[0])


ids = [f"T_c{c}_s{s}" for c in (0.65, 0.85, 1.05, 1.25, 1.5) for s in (9701, 9702)]
out = {}
for lam in ("0.5", "0.0625"):
    rows = {}
    for d in ids:
        an = git_json(f"qscan_anchor_theory_20261030/results/{d}__a1.0.json"); w = git_json(f"qscan_anchor_volve_20261029/results/dev/{d}.json")
        ha = curv(an["diag"]["J"][[str(l) for l in an["lams"]].index(lam)])
        hw = curv(w["diag"]["J"][[str(l) for l in w["lams"]].index(lam)])
        rows[d] = dict(H_a=ha, H_w=hw, ratio=(ha / hw) if ha and hw and hw > 0 else None)
    r = np.array([v["ratio"] for v in rows.values() if v["ratio"] is not None]); rng = np.random.default_rng(20261030)
    bs = [np.median(rng.choice(r, len(r))) for _ in range(2000)]
    out[lam] = dict(n=len(r), median=float(np.median(r)), ci=[float(np.percentile(bs, 5)), float(np.percentile(bs, 95))],
                    predicted=PRED, in_band=bool(4 <= np.median(r) <= 11), above2=bool(np.percentile(bs, 5) > 2), per_section=rows)
    print(f"lam={lam}: n={len(r)} median ratio={np.median(r):.2f} CI={np.round(out[lam]['ci'], 2)} "
          f"(pred {PRED}); ratios {np.round(np.sort(r), 1)}")
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "curvature.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
