"""Anchor experiment analysis (PLAN_ZH.md predictions 1-3) -> analysis/anchor.json.

Per weight: least squares v_sel = A u + K v_a + C over the 30 (dataset, anchor) runs (u = 1/c_true, v_a = 1/c_anchor,
v_sel = 1/c_sel); bootstrap over the 10 datasets (2000 resamples, seed 20261030, 90 % percentile) for K and A + K;
monotonicity of v_sel in v_a per dataset; agreement of the c_a = 1 run (lambda 0.5) with lib.qscan of experiment 12.
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]
FIT = json.loads((ROOT / "analysis" / "fit_existing.json").read_text(encoding="utf-8"))["dev (exp 12 sections)"]
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results").glob("*.json"))]
ids = sorted({r["id"] for r in recs}); lams = [str(l) for l in recs[0]["lams"]]


def fit(rs, lam):
    X = np.array([[1 / r["c_true"], 1 / r["c_anchor"], 1.0] for r in rs]); y = np.array([1 / r["qscan_a"][lam] for r in rs])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); return coef, float(np.sqrt(np.mean((y - X @ coef) ** 2)))


out = dict(n_runs=len(recs), datasets=ids, by_weight={})
rng = np.random.default_rng(20261030)
for lam in lams:
    (A, K, C), rms = fit(recs, lam); boot = []
    for _ in range(2000):
        pick = rng.choice(ids, len(ids)); rs = [r for d in pick for r in recs if r["id"] == d]
        if len({r["c_true"] for r in rs}) < 2:
            continue
        (a, k, _), _ = fit(rs, lam); boot.append((k, a + k))
    boot = np.array(boot)
    mono = {}
    for d in ids:
        v = [1 / r["qscan_a"][lam] for r in sorted((r for r in recs if r["id"] == d), key=lambda r: -r["c_anchor"])]
        mono[d] = dict(v_by_increasing_va=v, monotone=bool(np.all(np.diff(v) >= -1e-12)), strict=bool(np.all(np.diff(v) > 1e-12)))
    conv = float(np.mean([np.mean(r["diag"]["converged"][lams.index(lam)]) for r in recs]))
    w_key = "W" + lam
    out["by_weight"][lam] = dict(A=float(A), K=float(K), C=float(C), rms=rms, A_plus_K=float(A + K),
                                 K_ci=[float(np.percentile(boot[:, 0], 5)), float(np.percentile(boot[:, 0], 95))],
                                 A_plus_K_ci=[float(np.percentile(boot[:, 1], 5)), float(np.percentile(boot[:, 1], 95))],
                                 slope_qscan_w_dev=FIT[w_key]["A"], monotone_share=float(np.mean([m["monotone"] for m in mono.values()])),
                                 strict_share=float(np.mean([m["strict"] for m in mono.values()])), converged=conv, per_dataset=mono)
e12 = {d: json.loads(subprocess.run(["git", "show", f"HEAD:地震衰减补偿研究/baselines_volve_20261029/results/posthoc_qscan/{d}.json"],
                                    cwd=REPO, capture_output=True, check=True).stdout)["qscan"] for d in ids}
a1 = {r["id"]: r["qscan_a"]["0.5"] for r in recs if r["c_anchor"] == 1.0}
out["anchor1_vs_exp12"] = dict(same=int(sum(abs(a1[d] - e12[d]) < 1e-9 for d in ids)), n=len(ids),
                               pairs={d: (a1[d], e12[d]) for d in ids})
(ROOT / "analysis" / "anchor.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for lam, s in out["by_weight"].items():
    print(f"lam={lam}: A={s['A']:.2f} K={s['K']:.2f} {np.round(s['K_ci'], 2)} C={s['C']:+.2f} rms={s['rms']:.3f} "
          f"A+K={s['A_plus_K']:.2f} {np.round(s['A_plus_K_ci'], 2)} (QScan-W dev slope {s['slope_qscan_w_dev']:.2f}) "
          f"monotone {s['monotone_share']:.2f} strict {s['strict_share']:.2f} converged {s['converged']:.3f}")
print("c_a = 1, lambda 0.5 vs exp 12 lib.qscan:", out["anchor1_vs_exp12"]["same"], "/", len(ids))
