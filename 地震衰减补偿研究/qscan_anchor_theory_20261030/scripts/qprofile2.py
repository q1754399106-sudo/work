"""Second profile test (PLAN_PROFILE2_ZH.md): two new absorption profiles, the two-anchor kappa as a diagnostic.

usage: qprofile2.py    worker; jobs claimed with results10/<job>.claim, resumable.
Same generator, scan and selection as qprofile.py (after its addendum: early-stop scan, no QScan-W), with
profiles "D2" = [150, 150, 60, 60] (lossy below 2.04 s, milder than D; uniform kappa 0.489) and
"M" = [200, 45, 200, 200] (lossy between 0.72 and 2.04 s; uniform kappa 0.595); new seeds 9981-9984,
c 0.7 / 1.0 / 1.4, anchors c_a 0.7 / 1.0 / 1.4 (default wavelet window), lambda 1/2.
-> results10/<profile>_<dataset>__a<c_a>.json
"""
import os, json
import qprofile as qp
from env15 import np, gen, bm, ROOT

qp.PROFILES.clear(); qp.PROFILES.update({"D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0]})
SEEDS = (9981, 9982, 9983, 9984); OUT = ROOT / "results10"; bm.MAXITER = 20000


def jobs():
    return [(p, f"c{c}_s{s}", c, s, f"a{a}") for s in SEEDS for p in qp.PROFILES for c in qp.CS for a in (0.7, 1.0, 1.4)]


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
