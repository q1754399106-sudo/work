import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results"
R = []
for m in ("M1", "M10", "M11"):
    R += json.loads((OUT / f"diag_aligned_{m}.json").read_text(encoding="utf-8"))
print("repro max |eg - stored|", max(abs(r["eg_10-60"] - r["stored_eg"]) for r in R))
METH = ("method", "split", "direct"); SC = ("events_aligned", "thin_aligned", "dense_aligned"); FB = np.arange(10, 130, 10)


def cellmed(sel, k):
    """median of cell means (scene, M, snr), as in the step 2 tables"""
    cells = {}
    for r in sel:
        cells.setdefault((r["scene"], r["M"], r["snr"]), []).append(r[k])
    return float(np.median([np.mean(v) for v in cells.values()]))


for name, sel in [("all", R)] + [(s, [r for r in R if r["scene"] == s]) for s in SC]:
    print(f"\n=== {name}: E_G, then M_G (signed), shape, noise, noise-free bias; band 10-60")
    for m in METH:
        x = [r for r in sel if r["m"] == m]; f = lambda k: cellmed(x, k + "_10-60")
        print(f"{m:7s} eg {f('eg'):.3f} sgn_bias {f('sg'):+.3f} |bias| {f('mg'):.3f} shape {f('shape'):.3f} noise {f('noise'):.3f} bias_nf {f('mgs'):+.3f} ba {f('ba'):.3f} bg {f('bg'):.3f} lam {cellmed(x, 'lam'):.2e}" + (f" fc {cellmed(x, 'fc'):.1f} boundfrac {cellmed(x, 'bound_frac'):.3f} far {cellmed(x, 'bound_frac_far'):.3f}" if m == "method" else ""))
    print("band dependence: E_G and signed bias")
    for b in ("10-40", "10-60", "10-100", "none"):
        print(f"  {b:7s} " + " | ".join(f"{m} {cellmed([r for r in sel if r['m'] == m], 'eg_' + b):.3f} {cellmed([r for r in sel if r['m'] == m], 'mgs_' + b):+.3f}" for m in METH))
print("\nTransfer per band, noise-free, gB/gA, median over columns: f lo | method | split | direct")
for s in SC:
    print(s)
    for j, lo in enumerate(FB):
        row = [f"{lo:4d}"]
        for m in METH:
            x = [r for r in R if r["m"] == m and r["scene"] == s]; ga = np.median([r["ga"][j] for r in x]); gb = np.median([r["gb"][j] for r in x]); row.append(f"{ga:6.3f} {gb:6.3f} {gb / ga if abs(ga) > 1e-3 else float('nan'):6.2f}")
        print("  " + " | ".join(row))
for s in SC:
    for M in ("M1", "M10", "M11"):
        print(s, M, " | ".join(f"{m} eg {cellmed([r for r in R if r['m'] == m and r['scene'] == s and r['M'] == M], 'eg_10-60'):.3f} bias_nf {cellmed([r for r in R if r['m'] == m and r['scene'] == s and r['M'] == M], 'mgs_10-60'):+.3f}" for m in METH))
