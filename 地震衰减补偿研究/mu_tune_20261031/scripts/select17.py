"""Development choice of mu_relative (PLAN_ZH.md section 3) -> CHOICE.json and analysis/dev_summary.json.
Per route: the grid value with the highest mean S over the 24 development datasets (ties: the larger mu)."""
import json
from pathlib import Path
import numpy as np
from lib17 import GRID

ROOT = Path(__file__).resolve().parents[1]
R = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "dev").glob("*.json"))]
assert len(R) == 24, len(R)


def table(rows):
    return {route: {f"{mu:g}": float(np.mean([r["scores"][f"{route}_mu{mu:g}"]["S"] for r in rows])) for mu in GRID}
            for route in ("SR", "QT")}


def best(tab):
    v = max(tab.values()); return max(float(m) for m, s in tab.items() if s >= v - 1e-12)


allt = table(R); choice = {route: best(allt[route]) for route in allt}
out = dict(n=len(R), mean_S=allt, IQ70=float(np.mean([r["scores"]["IQ70"]["S"] for r in R])), choice=choice,
           by_batch={b: {route: best(t) for route, t in table([r for r in R if r["structure"] in ks]).items()}
                     for b, ks in (("first", "ABC"), ("second", "DEF"))},
           by_structure={k: {route: best(t) for route, t in table([r for r in R if r["structure"] == k]).items()} for k in "ABCDEF"},
           by_c0={str(c): {route: best(t) for route, t in table([r for r in R if r["c0"] == c]).items()} for c in (0.75, 1.05)},
           wins_vs_mu10={route: {f"{mu:g}": int(sum(r["scores"][f"{route}_mu{mu:g}"]["S"] > r["scores"][f"{route}_mu10"]["S"] for r in R))
                                 for mu in GRID} for route in ("SR", "QT")})
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "dev_summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
(ROOT / "CHOICE.json").write_bytes((json.dumps(choice) + "\n").encode("utf-8"))
print(json.dumps(out, indent=1))
