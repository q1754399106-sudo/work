"""Summary of the stage 0 results. Reads results/*.json; writes results/trend_summary.json. Descriptive only."""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"
T = json.loads((OUT / "trend.json").read_text(encoding="utf-8"))
P = json.loads((OUT / "precision.json").read_text(encoding="utf-8"))

base = {}
for did, v in P.items():
    var, c0, r = did[0], float(did.split("_c")[1].split("_")[0]), int(did[-1])
    base[(var, c0, r)] = v

rows = {}
for key, row in T.items():
    s = float(key[1:].split("_")[0]); c0 = float(key.split("_c")[1].split("_")[0]); r = int(key[-1])
    for var, v in row.items():
        b = base[(var, c0, r)]
        for kind in ("noisy", "clean"):
            rows.setdefault((s, var, kind), []).append(dict(
                sr=v[kind]["log_err_sr"], hsr=v[kind]["log_err_hsr"],
                dsr=v[kind]["log_err_sr"] - b[kind]["log_err_sr"]))

summary = {}
print(f"{'s':>8} {'var':>3} {'kind':>6} {'n':>2} {'SR err':>7} {'SR d vs s=0':>11} {'min':>6} {'max':>6} {'HSR sd':>7}")
for (s, var, kind) in sorted(rows):
    a = rows[(s, var, kind)]; sr = np.array([x["sr"] for x in a]); d = np.array([x["dsr"] for x in a]); h = np.array([x["hsr"] for x in a])
    summary[f"{s:+.4f}_{var}_{kind}"] = dict(n=len(a), sr_mean=float(sr.mean()), dsr_mean=float(d.mean()), dsr_min=float(d.min()),
                                             dsr_max=float(d.max()), hsr_sd=float(h.std(ddof=1)))
    print(f"{s:+8.4f} {var:>3} {kind:>6} {len(a):2d} {sr.mean():+7.3f} {d.mean():+11.3f} {d.min():+6.3f} {d.max():+6.3f} {h.std(ddof=1):7.2f}")

print("\nno trend (s = 0), same c0 and r subset")
for var in "ZY":
    for kind in ("noisy", "clean"):
        e = [base[(var, c, r)][kind]["log_err_sr"] for c in (0.75, 1.05) for r in (0, 1)]
        print(var, kind, f"SR err mean {np.mean(e):+.3f}, range {min(e):+.3f} to {max(e):+.3f}")
        summary[f"+0.0000_{var}_{kind}"] = dict(n=len(e), sr_mean=float(np.mean(e)))

print("\nslope of SR log error against s, by kind (all variants)")
for kind in ("noisy", "clean"):
    xs, ys = [], []
    for (s, var, k), a in rows.items():
        if k == kind:
            for x in a:
                xs.append(s); ys.append(x["dsr"])
    b = float(np.polyfit(xs, ys, 1)[0]); summary[f"slope_{kind}"] = b
    print(kind, f"d(log c_SR)/ds = {b:.1f} per (1/Hz/s); expected near 1/0.033 = 30 if the trend maps fully into 1/c")

(OUT / "trend_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
print("\nper-dataset no-trend SR and HSR, noisy")
for did, v in P.items():
    n = v["noisy"]; print(did, f"truth {n['truth']:.3f} SR {n['sr']:.3f} HSR {n['hsr']:.3f} edge {n['at_grid_edge']} windows {n['windows']}")
