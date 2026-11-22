"""Post-hoc, descriptive: how much of the negative raw-record Q can the source ghost (gun depth 6 m, not in the synthetic)
explain? Per row of results/field79raw.json: predicted y_g = -slope over 12 Hz..fmax of ln[|sin(2 pi f d cos_w / c)|_far /
|sin(...)|_near], cos_w of the far row (stored) and of the near reference (same window, |x| ~ 200 m, ray table).
Prints 1000/Q_eff with y and with y - y_g. -> results/field79raw_srcghost.json"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import field79 as F, field79raw as W  # noqa: E402
D = 6.0; C = 1500.0
d = json.loads((F.ROOT / "results" / "field79raw.json").read_text()); rows = d["rows"]; zw = d["water_depth"]
R = F.Rays(1.0, zw); pt = {round(tc, 1): W.ptable(R, tc) for tc in F.CENT}; xr = 200.0; out = []
for r in rows:
    f = F.F[(F.F >= F.BAND[0]) & (F.F <= r["fmax"])]; cw0, _ = W.ghost(pt[round(r["tc"], 1)](xr)[0], zw)
    g = np.log(np.abs(np.sin(2 * np.pi * f * D * r["cw"] / C))) - np.log(np.abs(np.sin(2 * np.pi * f * D * cw0 / C)))
    r["y_g"] = float(-np.polyfit(f, g, 1)[0]); out.append(r)
blocks = [[r for r in out if r["block"] == b] for b in range(F.NBLK)]
fit = lambda key: F.fit1([dict(r, y=r[key]) for r in out])
res = dict(invQ_eff=fit("y"), invQ_eff_minus_ghost=F.fit1([dict(r, y=r["y"] - r["y_g"]) for r in out]),
           mean_y=float(np.mean([r["y"] for r in out])), mean_y_ghost=float(np.mean([r["y_g"] for r in out])),
           cw_range=[float(min(r["cw"] for r in out)), float(max(r["cw"] for r in out))])
print(json.dumps(res, indent=1)); (F.ROOT / "results" / "field79raw_srcghost.json").write_text(json.dumps(res, indent=1))
