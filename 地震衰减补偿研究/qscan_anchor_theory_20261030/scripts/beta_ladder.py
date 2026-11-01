"""Exploratory (BETA_ZH.md): which field-like ingredient raises the offset beta of the QScan selection?

usage: beta_ladder.py [variant ...]    worker; jobs claimed with results12/<job>.claim, resumable.
Experiment-14 recipe (structure line 10227, noise lines 10149/10273/10355, curtain model Q), c 0.8 / 1.0 / 1.25,
seeds 9995-9998; fixed-anchor QScan (default wavelet window, c_a 0.7 / 1.0 / 1.4, lambda 1/2, scan and vertex as
qprofile.scan). beta = v(1.0) - (1 - kappa_hat) u - kappa_hat, kappa_hat from the 0.7 / 1.4 pair.
Variants:
  R0  block-bootstrap F-11 T2 log reflectivity (gen.signal), noise x 1
  R1  whole real well segments of 11 Volve wells (gen.signal with wells.all_segments), noise x 1
  R2  as R0, noise x 2
-> results12/<variant>_c<c>_s<seed>__a<c_a>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from qprofile import scan

sys.path.insert(0, str(ROOT.parent / "specslope_q_volve_20261028" / "scripts"))
OUT = ROOT / "results12"; bm.MAXITER = 20000
STRUCT, NOISE = 10227, [10149, 10273, 10355]
CS = (0.8, 1.0, 1.25); SEEDS = (9995, 9996, 9997, 9998); ANCH = (0.7, 1.0, 1.4)
VARIANTS = {"R0": dict(segs=False, mult=1.0), "R1": dict(segs=True, mult=1.0), "R2": dict(segs=False, mult=2.0)}


def jobs(vs):
    return [(v, c, s, a) for s in SEEDS for v in vs for c in CS for a in ANCH]


if __name__ == "__main__":
    vs = [a for a in sys.argv[1:] if a in VARIANTS] or list(VARIANTS)
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(STRUCT, NOISE); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    segs = None
    for v, c, s, a in jobs(vs):
        job = f"{v}_c{c}_s{s}__a{a}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        cfg = VARIANTS[v]
        if cfg["segs"] and segs is None:
            import wells
            segs = [x for w in wells.all_segments().values() for x in w]
        t0 = time.monotonic(); g = gen.signal(ctx, c, s, segs if cfg["segs"] else None)
        D = g["x"] + cfg["mult"] * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
        wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * a)); res = scan(D, cm, Q, pp, lambda cc: wv)
        rec = dict(variant=v, c_true=c, u=1 / c, seed=s, c_anchor=a, v_anchor=1 / a, noise_mult=cfg["mult"],
                   level=g["level"], seconds=time.monotonic() - t0, **res)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", f"{rec['seconds']:.0f}s", flush=True)
