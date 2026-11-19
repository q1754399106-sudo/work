"""Verdict of PLAN_GATE_A3_ZH.md. Cell = (model, fp, scene, c, snr), mean over 6 seeds; metric 'ga' (noisy per-event G/A)."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
R = [r for m in ("M28", "M29", "M30", "M31") for fp in (60, 40) for r in json.loads((ROOT / "results" / ("a3_%s_f%d.json" % (m, fp))).read_text())]


def cells(m, q="ga"):
    d = {}
    for r in R:
        if r["method"] == m:
            d.setdefault((r["model"], r["fp"], r["scene"], r["c"], r["snr"]), []).append(r[q])
    return {k: float(np.mean(v)) for k, v in d.items()}


S = cells("sat"); out = {"n_cells": len(S)}
for name, m in (("F1", "split"), ("F2", "band3")):
    B = cells(m); d = np.array([S[k] - B[k] for k in S]); q25, q75 = np.percentile(d, [25, 75])
    out[name] = dict(median_diff=float(np.median(d)), q25=float(q25), q75=float(q75), max_abs=float(np.abs(d).max()))
    out[name + "_pass"] = bool(abs(np.median(d)) <= 0.02 and q25 >= -0.03 and q75 <= 0.03)
B = cells("tik0"); r = np.array([B[k] / S[k] for k in S])
out["F3"] = dict(median_ratio=float(np.median(r)), frac_ge2=float(np.mean(r >= 2))); out["F3_pass"] = bool(np.median(r) >= 5 and np.mean(r >= 2) >= 0.9)
B = cells("direct"); r = np.array([B[k] / S[k] for k in S])
out["F4"] = dict(median_ratio=float(np.median(r))); out["F4_pass"] = bool(np.median(r) >= 1.2)
out["all_pass"] = all(v for k, v in out.items() if k.endswith("_pass"))
d = {}
for m in ("sat", "direct", "split", "band3", "tik0"):
    for q in ("ga", "wave", "eg"):
        C = cells(m, q); d[f"{m}_{q}_median"] = float(np.median(list(C.values())))
        for c in (0.5, 0.7, 1.0, 1.4, 2.0):
            d[f"{m}_{q}_c{c}"] = float(np.median([v for k, v in C.items() if k[3] == c]))
ch = [np.log(r["c_hat"] / r["c"]) for r in R if r["method"] == "sat" and r["scene"] == "events_aligned" and r["c_hat"] == r["c_hat"] and r["c_hat"] > 0]
d["ln_chat_rms"] = float(np.sqrt(np.mean(np.square(ch)))); d["ln_chat_mean"] = float(np.mean(ch))
d["chat_fail"] = int(sum(1 for r in R if r["method"] == "sat" and r["scene"] == "events_aligned" and not (r["c_hat"] == r["c_hat"] and r["c_hat"] > 0)))
out["descriptive"] = d
(ROOT / "results" / "a3_verdict.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "descriptive"}, indent=1))
