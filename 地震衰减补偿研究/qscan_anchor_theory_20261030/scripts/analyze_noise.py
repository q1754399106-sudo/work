"""Noise test analysis (PLAN_NOISE_ZH.md predictions 1-4) -> analysis/noise.json.

Per noise level: QScan-W fit v_sel = A_W u + C (6 datasets); fixed-anchor fit v_sel = A u + K v_a + C per wavelet
window (6 datasets x 2 anchors). Bootstrap over datasets within the level (2000 resamples, seed 20261030, 90 %).
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results3").glob("*.json"))]
MULTS = sorted({r["noise_mult"] for r in recs}); KAPPA = {"w1.4": 0.675, "w2.6": 1.14}


def fit_w(rs):
    u = np.array([1 / r["c_true"] for r in rs]); v = np.array([1 / r["c_sel"] for r in rs])
    return float(np.polyfit(u, v, 1)[0])


def fit_a(rs):
    X = np.array([[1 / r["c_true"], 1 / r["c_anchor"], 1.0] for r in rs]); y = np.array([1 / r["c_sel"] for r in rs])
    A, K, C = np.linalg.lstsq(X, y, rcond=None)[0]; return float(A), float(K), float(C)


def boot(rs, fn):
    ids = sorted({r["id"] for r in rs}); rng = np.random.default_rng(20261030); out = []
    for _ in range(2000):
        pick = rng.choice(ids, len(ids)); sub = [r for d in pick for r in rs if r["id"] == d]
        if len({r["c_true"] for r in sub}) > 1:
            out.append(fn(sub))
    return np.array(out)


ci = lambda x: [float(np.percentile(x, 5)), float(np.percentile(x, 95))]
out = dict(levels={})
for m in MULTS:
    rm = [r for r in recs if r["noise_mult"] == m]; lv = {}
    rw = [r for r in rm if r["scan"] == "qscanw"]; b = boot(rw, fit_w)
    lv["qscanw"] = dict(n=len(rw), A_W=fit_w(rw), A_W_ci=ci(b), converged=float(np.mean([np.mean(r["diag"]["converged"]) for r in rw])),
                        sel=[(r["c_true"], r["seed"], r["c_sel"]) for r in rw])
    for w in ("w1.4", "w2.6"):
        ra = [r for r in rm if r["scan"].startswith(w)]; A, K, C = fit_a(ra); b = boot(ra, lambda s: fit_a(s)[:2])
        lv[w] = dict(n=len(ra), A=A, K=K, C=C, A_plus_K=A + K, A_ci=ci(b[:, 0]), K_ci=ci(b[:, 1]), A_plus_K_ci=ci(b.sum(1)),
                     converged=float(np.mean([np.mean(r["diag"]["converged"]) for r in ra])),
                     sel=[(r["c_true"], r["seed"], r["c_anchor"], r["c_sel"]) for r in ra])
    out["levels"][str(m)] = lv
L = out["levels"]; s = [L[str(m)]["qscanw"]["A_W"] for m in MULTS]; dev = lambda m: abs(L[str(m)]["w2.6"]["A_plus_K"] - 1)
out["predictions"] = {
    "1_qscanw_slope_monotone": bool(np.all(np.diff(s) > 0)), "1_qscanw_lowest_noise_in_0.8_1.2": bool(0.8 <= s[0] <= 1.2),
    "2_deep_sum_lowest_noise_within_0.12": bool(dev(MULTS[0]) <= 0.12), "2_deep_sum_better_than_highest": bool(dev(MULTS[0]) < dev(MULTS[-1])),
    "3_K_deep_within_0.15": all(abs(L[str(m)]["w2.6"]["K"] - KAPPA["w2.6"]) <= 0.15 for m in MULTS),
    "3_K_default_within_0.15": all(abs(L[str(m)]["w1.4"]["K"] - KAPPA["w1.4"]) <= 0.15 for m in MULTS),
    "4_default_sum_within_0.2": all(abs(L[str(m)]["w1.4"]["A_plus_K"] - 1) <= 0.2 for m in MULTS)}
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "noise.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for m in MULTS:
    lv = L[str(m)]; q = lv["qscanw"]
    print(f"noise x{m}: QScan-W A_W={q['A_W']:.2f} {np.round(q['A_W_ci'], 2)} conv={q['converged']:.2f}")
    for w in ("w1.4", "w2.6"):
        a = lv[w]
        print(f"   {w}: A={a['A']:+.2f} K={a['K']:.2f} {np.round(a['K_ci'], 2)} A+K={a['A_plus_K']:.2f} {np.round(a['A_plus_K_ci'], 2)} conv={a['converged']:.2f}")
print(out["predictions"])
