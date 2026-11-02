"""Validity-boundary test (PLAN_BOUND_ZH.md): where does the uniform kappa formula hold?

usage: bound.py    worker; jobs claimed with results26/<job>.claim, resumable.
Same generator, scan and selection as qprofile2.py (early-stop scan, lambda 1/2, default wavelet window), seven new
absorption profiles (layer boundaries 0.72 / 2.04 / 2.82 s as the curtain model), new seeds 9991-9992,
c 0.7 / 1.0 / 1.4, anchors c_a 0.7 / 1.0 / 1.4.
-> results26/<profile>_<dataset>__a<c_a>.json
"""
import os, json
import qprofile as qp
from env15 import np, gen, bm, ROOT

PROFILES = {"R1": [100.0, 100.0, 100.0, 100.0], "R2": [60.0, 60.0, 150.0, 150.0], "R3": [150.0, 40.0, 150.0, 150.0],
            "R4": [50.0, 90.0, 250.0, 250.0], "X1": [120.0, 120.0, 70.0, 70.0], "X2": [300.0, 300.0, 60.0, 60.0],
            "X3": [150.0, 150.0, 150.0, 50.0]}
qp.PROFILES.clear(); qp.PROFILES.update(PROFILES)
SEEDS = (9991, 9992); OUT = ROOT / "results26"; bm.MAXITER = 20000


def jobs():
    return [(p, f"c{c}_s{s}", c, s, f"a{a}") for s in SEEDS for p in PROFILES for c in qp.CS for a in (0.7, 1.0, 1.4)]


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = gen.prepare(qp.STRUCT, qp.NOISE)
    for prof, did, c, s, sc in jobs():
        job = f"{prof}_{did}__{sc}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = qp.run(prof, did, c, s, sc, base)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", "edge" if rec["edge"] else "", f"{rec['seconds']:.0f}s", flush=True)
