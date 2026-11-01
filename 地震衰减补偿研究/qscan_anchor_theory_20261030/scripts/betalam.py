"""beta(lambda) (PLAN_BETALAM_ZH.md): selections at lambda 1/16, 1/8, 1/4 (+1/2) on Volve semi-synthetics and field lines.

usage: betalam.py    worker; jobs claimed with results23/<job>.claim, resumable.
Semi: base profile (gen recipe of qprofile.py, structure 10227, noise x 0.5), seeds 9971-9974 x c 0.7/1.0/1.4, anchor at
the truth (c_a = c, so v_sel = u + beta), lambda 1/16, 1/8, 1/4, 1/2. Field: the 26 SR inlines of field_kappa2.py, the
current pipeline (default wavelet window, c_a = 1), lambda 1/16, 1/8, 1/4 (1/2: results11). Scan: qprofile.scan with
qprofile.LAM set to lambda (lambda enters the solve and the L1 term of J).
-> results23/<job>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
import qprofile as qp
from field_kappa2 import LINES, L13

OUT = ROOT / "results23"; bm.MAXITER = 20000
SLAMS = (0.0625, 0.125, 0.25, 0.5); FLAMS = (0.0625, 0.125, 0.25)


def run_semi(c, s, lam, base):
    qp.LAM = lam; cm = dict(base["cm"]); Q = np.asarray(cm["Q"], float)
    g = gen.signal(base, c, s); D = g["x"] + qp.MULT * g["level"] * base["noise"][0]; pp = sc11.prepare_sc(D, cm)
    wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * c)); t0 = time.monotonic(); r = qp.scan(D, cm, Q, pp, lambda cc: wv)
    return dict(kind="semi", id=f"c{c}_s{s}", c_true=c, u=1 / c, c_anchor=c, lam=lam, seconds=time.monotonic() - t0, **r)


def run_field(ln, lam, cm, Q):
    qp.LAM = lam; D = L13.read_inline(ln); pp = sc11.prepare_sc(D, cm); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * 1.0))
    t0 = time.monotonic(); r = qp.scan(D, cm, Q, pp, lambda cc: wv)
    return dict(kind="field", inline=ln, c_anchor=1.0, lam=lam, seconds=time.monotonic() - t0, **r)


def jobs():
    J = [(f"S_c{c}_s{s}_lam{lam}", ("s", c, s, lam)) for lam in SLAMS for s in (9971, 9972, 9973, 9974) for c in qp.CS]
    J += [(f"F_IL_{ln}_lam{lam}", ("f", ln, lam)) for lam in FLAMS for ln in LINES]
    return J


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = None; cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for job, spec in jobs():
        path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        if spec[0] == "s":
            base = base or gen.prepare(qp.STRUCT, qp.NOISE); rec = run_semi(*spec[1:], base)
        else:
            rec = run_field(*spec[1:], cm, Q)
        path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", "edge" if rec["edge"] else "", f"{rec['seconds']:.0f}s", flush=True)
