"""Exploratory (post hoc, no frozen plan): why does kappa on the deep-lossy profile D exceed the uniform-weight value?

usage: dtrunc.py    worker; jobs claimed with results8/<job>.claim, resumable.
Same data as qprofile.py (profile D, seeds 9971-9972, c 0.7 / 1.0 / 1.4). Fixed-anchor QScan (c_a 0.7 / 1.4,
lambda 1/2) with the analysis window cut to 1.2-T1 s (T1 = 2.2, 2.6) as in trunc.py; wavelet from the full section.
If the deep part carries little Q information, K changes little between T1 = 2.6 and 3.0 (uniform weights predict
0.448 -> 0.294); if weights are uniform, K follows a0<a>/<a^2> of the cut window (2.2: 0.749, 2.6: 0.448).
-> results8/D_<dataset>__a<c_a>_t<T1>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from qprofile import context, scan, MULT, STRUCT, NOISE
from trunc import cut

OUT = ROOT / "results8"; bm.MAXITER = 20000


def jobs():
    return [(f"c{c}_s{s}", c, s, a, t) for s in (9971, 9972) for t in (2.2, 2.6) for c in (0.7, 1.0, 1.4) for a in (0.7, 1.4)]


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = gen.prepare(STRUCT, NOISE); ctx = context(base, "D"); cm = ctx["cm"]
    Q = np.asarray(cm["Q"], float)
    for did, c, s, ca, t1 in jobs():
        job = f"D_{did}__a{ca}_t{t1}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic(); g = gen.signal(ctx, c, s); D = g["x"] + MULT * g["level"] * ctx["noise"][0]
        pp = sc11.prepare_sc(D, cm); Dt, ct, pt = cut(D, cm, pp, t1); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * ca))
        res = scan(Dt, ct, Q, pt, lambda cc: wv)
        rec = dict(profile="D", id=did, c_true=c, u=1 / c, seed=s, c_anchor=ca, v_anchor=1 / ca, t1=t1,
                   seconds=time.monotonic() - t0, **res)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", f"{rec['seconds']:.0f}s", flush=True)
