"""usage: dev_sat_look.py [q]"""
import gzip, json, sys
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results"; quick = len(sys.argv) > 1; Ms = ("M1", "M10", "M11")
R = []; Base = []
for M in Ms:
    p = OUT / (f"dev_sat_{M}" + ("_q" if quick else "") + ".json")
    if p.exists():
        R += json.loads(p.read_text(encoding="utf-8"))
    Base += json.loads(gzip.open(OUT / f"step2_{M}.json.gz", "rt", encoding="utf-8").read())
qs = {r["qf"] for r in R}; sc = {r["scene"] for r in R}; A = [r for r in Base if r["qf"] in qs and r["scene"] in sc] + R
vars_ = ["method", "split", "direct"] + sorted({r["m"] for r in R})


def tab(sel, title):
    print(f"\n{title}")
    idx = {}
    for r in sel:
        idx.setdefault((r["m"], r["scene"], r["M"], r["snr"], r["qf"]), []).append(r)
    for m in vars_:
        cells = [v for k, v in idx.items() if k[0] == m]
        if not cells:
            continue
        f = lambda k: float(np.median([np.mean([x[k] for x in v]) for v in cells]))
        print(f"  {m:9s} eg {f('eg'):.3f} cg {f('cg'):.3f} wave {f('wave'):.3f} mg {f('mg'):.3f} ma {f('ma'):.3f} lam {f('lam'):.2e} n={len(cells)}")


tab([r for r in A if r["qf"] == 1.0 and "aligned" in r["scene"]], "aligned Q correct")
for s in ("events", "thin", "dense"):
    tab([r for r in A if r["qf"] == 1.0 and r["scene"] == s + "_aligned"], f"{s} aligned Q correct")
if not quick:
    tab([r for r in A if r["qf"] != 1.0 and "aligned" in r["scene"]], "aligned Q wrong")
    tab([r for r in A if r["qf"] == 1.0 and "nmo" in r["scene"]], "NMO Q correct")
    tab([r for r in A if r["qf"] != 1.0 and "nmo" in r["scene"]], "NMO Q wrong")
    for s in ("events", "thin", "dense"):
        tab([r for r in A if r["scene"] == s + "_nmo"], f"{s} nmo all Q")
    for M in Ms:
        tab([r for r in A if r["qf"] == 0.8 and "aligned" in r["scene"] and r["M"] == M], f"{M} aligned Q x0.8")
print("\ntemplate bias events aligned Q correct (variants only):")
for m in vars_:
    x = [r["tga"] for r in R if r["m"] == m and r["scene"] == "events_aligned" and r["qf"] == 1.0]
    if x:
        print(f"  {m:9s} {np.mean(x):.3f}")
print("fallbacks:", {m: sum(r.get("fallback", 0) for r in R if r["m"] == m) for m in vars_})
