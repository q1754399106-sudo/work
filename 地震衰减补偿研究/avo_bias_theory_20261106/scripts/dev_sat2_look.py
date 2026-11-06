import gzip, json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results"; R = []
for M in ("M1", "M10", "M11"):
    p = OUT / f"dev_sat2_{M}.json.gz"
    if p.exists():
        R += json.loads(gzip.open(p, "rt", encoding="utf-8").read())
V = ("split", "method", "direct", "sat6", "sat10", "sat15", "sat25"); K = ("gacb", "gab", "gacf", "gaf", "egb", "eglb", "mgb", "mglcb", "cgb", "egf", "wave", "lam")
print("rows", len(R), "fallback", {v: sum(r["fallback"] for r in R if r["m"] == v) for v in V})


def tab(sel, title):
    idx = {}
    for r in sel:
        idx.setdefault((r["m"], r["scene"], r["M"], r["snr"], r["qf"]), []).append(r)
    print(f"\n{title}\n  {'':7s}" + "".join(f"{k:>8s}" for k in K) + "   n")
    for m in V:
        cells = [v for k, v in idx.items() if k[0] == m]
        if cells:
            print(f"  {m:7s}" + "".join((f"{np.median([np.mean([x[k] for x in v]) for v in cells]):8.3f}" if k != "lam" else f"{np.median([np.mean([x[k] for x in v]) for v in cells]):8.1e}") for k in K) + f" {len(cells):3d}")


al = [r for r in R if "aligned" in r["scene"]]; nm = [r for r in R if "nmo" in r["scene"]]; s = lambda X, c: [r for r in X if c(r)]
tab(s(al, lambda r: r["qf"] == 1), "aligned Q correct")
for sc in ("events", "thin", "dense"):
    tab(s(al, lambda r: r["qf"] == 1 and r["scene"].startswith(sc)), f"aligned Q correct {sc}")
tab(s(al, lambda r: r["qf"] != 1), "aligned Q wrong"); tab(s(al, lambda r: r["qf"] == .8), "aligned Q x0.8"); tab(s(al, lambda r: r["qf"] == 1.2), "aligned Q x1.2")
tab(nm, "NMO all"); tab(s(nm, lambda r: r["qf"] == 1), "NMO Q correct"); tab(s(nm, lambda r: r["qf"] != 1), "NMO Q wrong")
for sc in ("events", "thin", "dense"):
    tab(s(nm, lambda r: r["scene"].startswith(sc)), f"NMO {sc}")
for M in ("M1", "M10", "M11"):
    tab(s(al, lambda r: r["M"] == M and r["qf"] == 1), f"{M} aligned Q correct")
print("\nsigned noise-free per-event G/A error, band 10-60, aligned events Q correct")
for m in V:
    x = s(al, lambda r: r["m"] == m and r["qf"] == 1 and r["scene"].startswith("events"))
    if x:
        print(f"  {m:7s}" + " ".join(f"{np.mean([r[f'gac{j}b'] for r in x]):+.4f}" for j in range(3)) + "  full band " + " ".join(f"{np.mean([r[f'gac{j}f'] for r in x]):+.4f}" for j in range(3)))
