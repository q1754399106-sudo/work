"""Analysis-window test analysis (PLAN_TRUNC_ZH.md predictions 1-3) -> analysis/trunc.json.

QScan-W per T1: v_sel = A_W u + C (12 datasets); fixed anchor at T1 = 2.2 s: v_sel = A u + K v_a + C (24 scans).
Bootstrap over datasets (2000, seed 20261030, 90 %); A_W(3.0) - A_W(2.2) paired over datasets.
Reliability: share of T1 = 2.2 QScan-W picks at a grid edge or with an unconverged point among the pick and its neighbours.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results4").glob("*.json"))]
T1S = [3.0, 2.6, 2.2]; KAPPA = 0.804; IDS = sorted({r["id"] for r in recs})
W = {t: {r["id"]: r for r in recs if r["scan"] == f"w_t{t}"} for t in T1S}
AN = [r for r in recs if r["scan"].startswith("a")]


def fit_w(ids, t):
    u = np.array([1 / W[t][i]["c_true"] for i in ids]); v = np.array([1 / W[t][i]["c_sel"] for i in ids])
    return float(np.polyfit(u, v, 1)[0]), float(np.polyfit(u, v, 1)[1])


def fit_a(ids):
    rs = [r for d in ids for r in AN if r["id"] == d]
    X = np.array([[1 / r["c_true"], 1 / r["c_anchor"], 1.0] for r in rs]); y = np.array([1 / r["c_sel"] for r in rs])
    return np.linalg.lstsq(X, y, rcond=None)[0]


def boot(fn):
    rng = np.random.default_rng(20261030); out = []
    for _ in range(2000):
        ids = list(rng.choice(IDS, len(IDS)))
        if len({W[3.0][i]["c_true"] for i in ids}) > 1:
            out.append(fn(ids))
    return np.array(out)


def unreliable(r):
    j = int(np.argmin(r["diag"]["J"])); n = len(r["diag"]["J"]); conv = r["diag"]["converged"]
    return j in (0, n - 1) or not all(conv[k] for k in range(max(j - 1, 0), min(j + 2, n)))


ci = lambda x: [float(np.percentile(x, 5)), float(np.percentile(x, 95))]
out = dict(n_datasets=len(IDS), qscanw={}, anchored={})
for t in T1S:
    A, C = fit_w(IDS, t); b = boot(lambda ids: fit_w(ids, t)[0])
    out["qscanw"][str(t)] = dict(A_W=A, C=C, A_W_ci=ci(b), unreliable=float(np.mean([unreliable(r) for r in W[t].values()])),
                                 converged=float(np.mean([np.mean(r["diag"]["converged"]) for r in W[t].values()])),
                                 sel=sorted((r["c_true"], r["seed"], r["c_sel"]) for r in W[t].values()))
d = boot(lambda ids: fit_w(ids, 3.0)[0] - fit_w(ids, 2.2)[0]); out["qscanw"]["diff_3.0_minus_2.2"] = dict(value=out["qscanw"]["3.0"]["A_W"] - out["qscanw"]["2.2"]["A_W"], ci=ci(d))
A, K, C = fit_a(IDS); b = boot(fit_a)
out["anchored"] = dict(A=float(A), K=float(K), C=float(C), A_plus_K=float(A + K), A_ci=ci(b[:, 0]), K_ci=ci(b[:, 1]), A_plus_K_ci=ci(b[:, 0] + b[:, 1]),
                       converged=float(np.mean([np.mean(r["diag"]["converged"]) for r in AN])),
                       sel=sorted((r["c_true"], r["seed"], r["c_anchor"], r["c_sel"]) for r in AN))
q = out["qscanw"]; a22 = q["2.2"]["A_W"]
out["predictions"] = {"1_control_A_W_3.0_in_1.1_1.5": bool(1.1 <= q["3.0"]["A_W"] <= 1.5),
                      "2_deep": ("undetermined (reliability)" if q["2.2"]["unreliable"] > 0.25 else
                                 "supported" if a22 <= 1.12 and q["diff_3.0_minus_2.2"]["value"] >= 0.15 else
                                 "rejected" if a22 >= 1.2 else "undetermined"),
                      "3_K_within_0.804_pm_0.10": bool(abs(K - KAPPA) <= 0.10), "3_A_plus_K_within_1_pm_0.2": bool(abs(A + K - 1) <= 0.2)}
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "trunc.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for t in T1S:
    s = q[str(t)]; print(f"T1 {t}: A_W={s['A_W']:.2f} {np.round(s['A_W_ci'], 2)} C={s['C']:+.2f} unreliable={s['unreliable']:.2f} conv={s['converged']:.2f}")
print("diff 3.0-2.2:", round(q["diff_3.0_minus_2.2"]["value"], 2), np.round(q["diff_3.0_minus_2.2"]["ci"], 2))
a = out["anchored"]; print(f"anchored T1 2.2: A={a['A']:+.2f} K={a['K']:.2f} {np.round(a['K_ci'], 2)} A+K={a['A_plus_K']:.2f} {np.round(a['A_plus_K_ci'], 2)}")
print(out["predictions"])
