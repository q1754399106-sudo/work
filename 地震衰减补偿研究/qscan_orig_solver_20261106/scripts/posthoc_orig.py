"""Post-hoc checks (not in PLAN_ZH.md): no new scans, reads ../results only -> ../analysis/posthoc.json."""
import json, glob, sys
from pathlib import Path
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
from analyze_orig import load, slope, kappa, num, boot

T = load(); qts = sorted({k[1] for k in T}); seeds = sorted({k[2] for k in T})
TS = {"W1": 0.1, "W2": 0.5, "W3": 1.0}; out = {}
grid = np.exp(np.linspace(np.log(1 / 400), np.log(1 / 20), 21))
out["grid_step_K_units"] = {str(q): float(np.interp(1 / q, grid[:-1], np.diff(grid)) / 0.01) for q in qts}

# 1. sum rule A + K (never scored by the plan)
for crit in ("v1", "v2"):
    for w in TS:
        fn = lambda ss: slope(T, f"{w}_0.00", crit, qts, ss) + kappa(T, w, crit, qts, ss)
        out[f"{crit}:{w}:A+K"] = dict(val=fn(seeds), ci=boot(fn, seeds))
    # 2. slope relative to the known wavelet
    for w in TS:
        fn = lambda ss: slope(T, f"{w}_0.00", crit, qts, ss) - slope(T, "known", crit, qts, ss)
        out[f"{crit}:{w}:A-Aknown"] = dict(val=fn(seeds), ci=boot(fn, seeds))
    # 3. mean bias of v_sel from u, per Q_true
    for key in ["known"] + [f"{w}_0.00" for w in TS]:
        out[f"{crit}:{key}:v_bias"] = {str(q): float(np.nanmedian([num(T[(key, q, s)][crit]) - 1 / q for s in seeds])) for q in qts}

# 4. K per Q_true; spread of per-pair K
for crit in ("v1", "v2"):
    for w in TS:
        out[f"{crit}:{w}:K_by_Q"] = {}
        for q in qts:
            d = [(num(T[(f"{w}_0.01", q, s)][crit]) - num(T[(f"{w}_0.00", q, s)][crit])) / 0.01 for s in seeds]
            out[f"{crit}:{w}:K_by_Q"][str(q)] = float(np.nanmedian(d))

# 5. measured peak frequency vs the Zhang-Ulrych prediction from f0 = 30 Hz, and the v_a it implies
for w, ts in TS.items():
    rows = {}
    for q in qts:
        fp = float(np.median([num(T[(f"{w}_0.00", q, s)]["fp"]) for s in seeds]))
        fw = 30.0; va_impl = (2 / fp - 2 * fp / fw**2) / (np.pi * ts)
        f = 30.0
        for _ in range(200):
            f = 1 / np.sqrt(1 / 30.0**2 + np.pi * ts * (1 / q) / (2 * f))
        rows[str(q)] = dict(fp=fp, fp_theory=float(f), v_true=1 / q, v_implied=float(va_impl))
    out[f"fp:{w}"] = rows
(ROOT / "analysis" / "posthoc.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for k, v in out.items():
    print(k, json.dumps(v, ensure_ascii=False))
