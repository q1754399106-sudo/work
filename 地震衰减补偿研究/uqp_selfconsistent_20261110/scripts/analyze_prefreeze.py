"""Summary of pilot_prefreeze.py (scenario mean over the 9 Q multipliers; units: dataset / FD model)."""
import json, glob
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results"
mean = lambda recs, k: float(np.mean([r[k] for r in recs]))
print("P2 semi-synthetic (24 datasets): scenario means and UQMW_K9 - STK")
for v in ("cal", "nd"):
    rows = []
    for p in sorted(glob.glob(str(R / "pilot_prefreeze" / f"semi_*_{v}.json"))):
        a = json.loads(Path(p).read_text(encoding="utf-8")); b = json.loads((R / "pilot_semi" / f"{a['id']}_{v}.json").read_text(encoding="utf-8"))["records"]
        rows.append((mean(b, "UQW"), mean(a["records"], "STK"), mean(b, "UQMW_K9")))
    rows = np.array(rows); d = rows[:, 2] - rows[:, 1]; rng = np.random.default_rng(20261116)
    bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
    print(f"  {v}: n {len(rows)}  UQW {rows[:,0].mean():.3f}  STK {rows[:,1].mean():.3f}  UQMW_K9 {rows[:,2].mean():.3f}  "
          f"UQMW-STK {d.mean():+.3f} [{np.quantile(bs,.05):+.3f}, {np.quantile(bs,.95):+.3f}] {int((d>0).sum())}/{len(d)}  STK-UQW {np.mean(rows[:,1]-rows[:,0]):+.3f}")
print("P1 FD H1-H6 (model = mean of 3 realizations)")
KEYS = ("UQW", "STK", "UQMW_K9_flat", "UQMW_K9_steer", "MW0_K9_flat", "MW0_K9_steer")
for lv in ("14dB", "30dB"):
    T = {}
    for p in glob.glob(str(R / "pilot_prefreeze" / f"fd_*_{lv}_*.json")):
        a = json.loads(Path(p).read_text(encoding="utf-8")); T.setdefault(a["model"], []).append([mean(a["records"], k) for k in KEYS] + [a["records"][4]["MW0_K9_steer"], a["records"][4]["MW0_K9_flat"], a["dip_abs_p95_ms"]])
    A = np.array([np.mean(T[h], axis=0) for h in sorted(T)])
    print(f"  {lv}: " + "  ".join(f"{k} {A[:,i].mean():.3f}" for i, k in enumerate(KEYS)))
    for name, i, j in (("steer-flat (UQMW)", 3, 2), ("UQMW_steer-UQW", 3, 0), ("UQMW_steer-STK", 3, 1), ("STK-UQW", 1, 0)):
        d = A[:, i] - A[:, j]; print(f"      {name}: {d.mean():+.3f}  {int((d>0).sum())}/{len(d)}  per model {np.round(d,3).tolist()}")
    print(f"      m=1 MW0 steer {A[:,6].mean():.3f} flat {A[:,7].mean():.3f}; dip |p| 95% {A[:,8].mean():.2f} ms/trace")
