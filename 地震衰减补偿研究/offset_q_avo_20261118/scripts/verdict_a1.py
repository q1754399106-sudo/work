"""Verdict of PLAN_GATE_A1_ZH.md section 3, and descriptive tables."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; GM = ["M16", "M17", "M18", "M19"]
R = [r for m in GM for r in json.loads((ROOT / "results" / ("a1_%s.json" % m)).read_text())]


def lerr(rows):
    return np.array([np.log(r["c_hat"] / r["c"]) if (r["c_hat"] == r["c_hat"] and r["c_hat"] > 0) else np.nan for r in rows])


def sel(kind="avo", snr=24, v=1.0, n=50, **k):
    return [r for r in R if r["kind"] == kind and r["snr"] == snr and r["v"] == v and r["n"] == n and all(r[a] == b for a, b in k.items())]


def summ(x):
    ok = x[np.isfinite(x)]
    return dict(k=len(x), nan=int(np.isnan(x).sum()), mean=float(ok.mean()) if len(ok) else None, sd=float(ok.std()) if len(ok) else None,
                rms=float(np.sqrt(np.mean(ok ** 2))) if len(ok) else None)


out = {}; x = lerr(sel()); s = summ(x)
out["A1_rms"] = s["rms"]; out["A1_pass"] = bool(s["rms"] <= 0.10); out["main"] = s
cells = {}
for m in GM:
    for c in (0.7, 1.0, 1.4):
        for fp in (40.0, 60.0):
            cells["%s_c%.1f_f%g" % (m, c, fp)] = summ(lerr(sel(model=m, c=c, fp=fp)))
worst = max(cells.items(), key=lambda kv: abs(kv[1]["mean"]) if kv[1]["mean"] is not None else 9)
out["A2_worst_cell"] = worst[0]; out["A2_worst_absmean"] = abs(worst[1]["mean"]); out["A2_pass"] = bool(out["A2_worst_absmean"] <= 0.10)
out["A3_bad_frac"] = float(np.mean(~(np.isfinite(x)))); out["A3_pass"] = bool(out["A3_bad_frac"] <= 0.02)
out["cells"] = cells
d = {}
for kind, snr, v in [("avo", 24, 1.0), ("avo", 30, 1.0), ("avo", 18, 1.0), ("flat", 24, 1.0), ("avo", 24, 1.02), ("avo", 24, 0.98)]:
    for n in (10, 50, 200):
        d["%s_%d_v%.2f_N%d" % (kind, snr, v, n)] = summ(lerr(sel(kind, snr, v, n)))
        for fp in (40.0, 60.0):
            d["%s_%d_v%.2f_N%d_f%g" % (kind, snr, v, n, fp)] = summ(lerr(sel(kind, snr, v, n, fp=fp)))
for c in (0.7, 1.0, 1.4):
    d["by_c_%.1f" % c] = summ(lerr(sel(c=c)))
for m in GM:
    d["by_model_" + m] = summ(lerr(sel(model=m)))
d["nwin_median"] = float(np.median([r["nwin"] for r in sel()]))
out["descriptive"] = d
(ROOT / "results" / "a1_verdict.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: out[k] for k in out if k not in ("cells", "descriptive")}, indent=1))
