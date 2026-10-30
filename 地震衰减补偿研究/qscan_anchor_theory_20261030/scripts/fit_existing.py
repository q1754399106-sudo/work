"""Post-hoc check of the linear model in 1/c space on existing results (experiments 12 and 14).

Theory (THEORY_ZH.md): with u = 1/c_true, v = 1/c_sel, v_a = 1/c_anchor,
    v_sel = (1 - kappa) u + kappa v_a + beta
anchored QScan (v_a = 1): slope 1 - kappa < 1, intercept kappa + beta; QScan-W (wavelet re-estimated per c): kappa = 0,
slope 1, intercept beta(lambda). Fits v_sel = A u + C by least squares per method and data group; the log-space slope
implied by the fit is d log c_sel / d log c_true = A u / v.
Reads results from git HEAD (the experiment directories are archived).
-> analysis/fit_existing.json
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]; BASE = "地震衰减补偿研究"


def git_json(path):
    return json.loads(subprocess.run(["git", "show", f"HEAD:{BASE}/{path}"], cwd=REPO, capture_output=True, check=True).stdout)


def ls(path):
    out = subprocess.run(["git", "ls-tree", "--name-only", f"HEAD:{BASE}/{path}"], cwd=REPO, capture_output=True, check=True)
    return [n for n in out.stdout.decode().split() if n.endswith(".json")]


E14 = "qscan_anchor_volve_20261029/results"; E12 = "baselines_volve_20261029/results/posthoc_qscan"
dev = [git_json(f"{E14}/dev/{n}") for n in ls(f"{E14}/dev")]
e12 = {r["id"]: r for r in (git_json(f"{E12}/{n}") for n in ls(E12)) if r["id"] in {d["id"] for d in dev}}
test = [git_json(f"{E14}/test/{n}") for n in ls(f"{E14}/test")]
groups = {"dev (exp 12 sections)": [dict(c_true=d["c_true"], qscan=e12[d["id"]]["qscan"], sr=e12[d["id"]]["c_sr"], **{f"W{k}": v for k, v in d["qscan_w"].items()}) for d in dev],
          "T": [dict(c_true=r["c_true"], qscan=r["qscan"], sr=r["c_sr"], **{f"W{k}": v for k, v in r["qscan_w"].items()}) for r in test if r["group"] == "T"],
          "S": [dict(c_true=r["c_true"], qscan=r["qscan"], sr=r["c_sr"], **{f"W{k}": v for k, v in r["qscan_w"].items()}) for r in test if r["group"] != "T"]}
methods = ["qscan", "W0.5", "W0.25", "W0.0625", "W0.015625", "sr"]
out = {}
for g, rs in groups.items():
    u = np.array([1 / r["c_true"] for r in rs]); out[g] = {}
    for m in methods:
        v = np.array([1 / r[m] for r in rs]); A, C = np.polyfit(u, v, 1); res = v - (A * u + C)
        edge = float(np.mean([r[m] <= 0.4 + 1e-9 or r[m] >= 2.5 - 1e-9 for r in rs]))
        logslope = float(np.polyfit(np.log(1 / u), np.log(1 / v), 1)[0])
        out[g][m] = dict(n=len(u), A=float(A), C=float(C), rms_res=float(np.sqrt(np.mean(res ** 2))), mean_shift=float(np.mean(v - u)),
                         log_slope=logslope, log_slope_from_fit=float(np.mean(A * u / (A * u + C))), at_grid_edge=edge)
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "fit_existing.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for g, d in out.items():
    print(g)
    for m, s in d.items():
        print(f"  {m:<10} A={s['A']:5.2f} C={s['C']:+5.2f} rms={s['rms_res']:.3f} mean(v-u)={s['mean_shift']:+.3f} "
              f"logslope={s['log_slope']:.2f} (from fit {s['log_slope_from_fit']:.2f}) edge={s['at_grid_edge']:.2f}")
