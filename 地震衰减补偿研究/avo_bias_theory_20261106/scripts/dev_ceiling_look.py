import json
from pathlib import Path
import numpy as np
OUT = Path(__file__).resolve().parents[1] / "results"; R = []
for M in ("M1", "M10", "M11"):
    R += json.loads((OUT / f"dev_ceiling_{M}.json").read_text(encoding="utf-8"))
K = ("eg", "cg", "wave", "ga", "egl"); cells = sorted({(r["scene"], r["qf"]) for r in R})
for sc, qf in cells:
    print(f"\n{sc} Q x{qf}: mean over models, columns")
    x = [r for r in R if r["scene"] == sc and r["qf"] == qf]
    for m in ("sat", "split"):
        lams = sorted({r["lam"] for r in x if r["m"] == m}, reverse=True); tab = {}
        for lam in lams:
            tab[lam] = {k: float(np.mean([np.mean(r[k]) for r in x if r["m"] == m and r["lam"] == lam])) for k in K}
        b = min(lams, key=lambda l: tab[l]["eg"]); bw = min(lams, key=lambda l: tab[l]["wave"])
        print(f"  {m:6s} oracle-E_G lam {b:.1e}: " + " ".join(f"{k} {tab[b][k]:.3f}" for k in K) + f"  | oracle-wave lam {bw:.1e}: wave {tab[bw]['wave']:.3f} eg {tab[bw]['eg']:.3f}")
        if m == "sat":
            ch = []
            for M in ("M1", "M10", "M11"):
                y = sorted([r for r in x if r["m"] == "sat" and r["M"] == M], key=lambda r: -r["lam"])
                for c in range(len(y[0]["eg"])):
                    pick = next((r for r in y if r["res"][c] <= 1), y[-1]); ch.append({k: pick[k][c] for k in K} | {"lam": pick["lam"]})
            print(f"  {'discr':6s} median lam {np.median([z['lam'] for z in ch]):.1e}: " + " ".join(f"{k} {np.mean([z[k] for z in ch]):.3f}" for k in K))
    y = [r for r in x if r["m"] == "band3"]; print(f"  {'band3':6s} " + " ".join(f"{k} {np.mean([np.mean(r[k]) for r in y]):.3f}" for k in K) + "  max wave " + f"{max(max(r['wave']) for r in y):.2f}")
