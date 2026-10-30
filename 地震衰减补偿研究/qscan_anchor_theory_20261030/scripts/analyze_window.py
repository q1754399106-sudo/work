"""Anchor-window analysis (PLAN_WINDOW_ZH.md criteria 1-4) -> analysis/window.json.

Per wavelet window: least squares v_sel = A u + K v_a + C over the 20 (dataset, anchor) scans; bootstrap over the
10 datasets (2000 resamples, seed 20261030, 90 % percentile) for A, K, A + K. The 1.2-1.6 s window is taken from the
anchor experiment (lambda 0.5, all three anchors) for criterion 2.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRED = {"1.4": 0.675, "2.0": 0.96, "2.6": 1.14}
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results2").glob("*.json"))]
base = json.loads((ROOT / "analysis" / "anchor.json").read_text(encoding="utf-8"))["by_weight"]["0.5"]
ids = sorted({r["id"] for r in recs})


def fit(rs):
    X = np.array([[1 / r["c_true"], 1 / r["c_anchor"], 1.0] for r in rs]); y = np.array([1 / r["c_sel"] for r in rs])
    return np.linalg.lstsq(X, y, rcond=None)[0]


out = dict(baseline_1p4=dict(A=base["A"], K=base["K"], A_plus_K=base["A_plus_K"], K_ci=base["K_ci"]), windows={})
for tc in ("2.0", "2.6"):
    rs = [r for r in recs if str(r["window"][2]) == tc]; A, K, C = fit(rs); rng = np.random.default_rng(20261030); bs = []
    for _ in range(2000):
        pick = rng.choice(ids, len(ids)); sub = [r for d in pick for r in rs if r["id"] == d]
        if len({r["c_true"] for r in sub}) > 1:
            bs.append(fit(sub)[:2])
    bs = np.array(bs); ci = lambda x: [float(np.percentile(x, 5)), float(np.percentile(x, 95))]
    out["windows"][tc] = dict(n=len(rs), A=float(A), K=float(K), C=float(C), A_plus_K=float(A + K), A_ci=ci(bs[:, 0]),
                              K_ci=ci(bs[:, 1]), A_plus_K_ci=ci(bs.sum(1)), kappa_pred=PRED[tc],
                              converged=float(np.mean([np.mean(r["diag"]["converged"]) for r in rs])),
                              at_edge=float(np.mean([r["c_sel"] <= 0.4 + 1e-9 or r["c_sel"] >= 2.5 - 1e-9 for r in rs])),
                              sel=[(r["id"], r["c_anchor"], r["c_sel"]) for r in rs])
w = out["windows"]
out["criteria"] = {"1_K_within_0.15": all(abs(w[t]["K"] - PRED[t]) <= 0.15 for t in w),
                   "2_K_increasing": bool(base["K"] < w["2.0"]["K"] < w["2.6"]["K"]),
                   "3_A_negative_2.6": bool(w["2.6"]["A"] < 0),
                   "4_A_plus_K_within_0.2": all(abs(w[t]["A_plus_K"] - 1) <= 0.2 for t in w)}
(ROOT / "analysis" / "window.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("baseline 1.4:", {k: round(v, 2) if isinstance(v, float) else v for k, v in out["baseline_1p4"].items()})
for t, s in w.items():
    print(f"window tc={t}: A={s['A']:+.2f} {np.round(s['A_ci'], 2)} K={s['K']:.2f} {np.round(s['K_ci'], 2)} (pred {s['kappa_pred']}) "
          f"A+K={s['A_plus_K']:.2f} {np.round(s['A_plus_K_ci'], 2)} C={s['C']:+.2f} converged={s['converged']:.3f}")
print(out["criteria"])
