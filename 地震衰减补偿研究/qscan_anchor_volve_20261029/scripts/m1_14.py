"""Experiment 14: M1 of the SC-a outputs in the gather figures (descriptive) -> analysis/gathers_m1.json.

M1 = sc11.m1 (median trace correlation with the ideal target, 2.30-2.80 s), read from cache/gathers/<id>.npz written
by gathers14.py; per group, the mean M1 for each Q factor choice (true c, SR, QScan, QScan-W(2^-4), its calibrated form).
"""
import json
from env14 import np, sc11, ROOT

out = {}
for z in sorted((ROOT / "cache" / "gathers").glob("*.npz")):
    a = dict(np.load(z)); D, ideal = a.pop("Input"), a.pop("ideal"); grp = "T" if z.stem.startswith("T_") else "S"
    row = {"输入": float(sc11.m1(D, ideal)["M1"])} | {k: float(sc11.m1(v, ideal)["M1"]) for k, v in a.items()}
    out.setdefault(grp, {})[z.stem] = row
summ = {g: {k: float(np.mean([r[k] for r in rs.values()])) for k in next(iter(rs.values()))} for g, rs in out.items()}
(ROOT / "analysis" / "gathers_m1.json").write_text(json.dumps(dict(mean=summ, per_dataset=out), indent=1, ensure_ascii=False), encoding="utf-8")
for g, s in summ.items():
    print(g, len(out[g]), {k: round(v, 3) for k, v in s.items()})
