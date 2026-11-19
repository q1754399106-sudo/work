"""Anchor-point lasso solution per dataset (PLAN_ZH.md): energy envelope and non-zero density of r0(t).

usage: active_curv.py [--smoke]    worker; jobs claimed with results/<job>.claim, resumable.
One lasso solve per dataset at the current-pipeline anchor c_a = 1: wavelet lib.wavelets(D, TW, cm[Q*1]), dictionary at
v = 1 (model Q), lambda 1/2, cold start, FISTA limit 20000. Stored per atom time: e(t) = mean over live traces of X^2,
n(t) = share of live traces with |X| > 1e-8, plus t, a(t) (model Q), a0 = a(1.4 s). Only the recorded section D, the
model Q and the anchor wavelet enter; the generator is called with the true c only to rebuild the archived data.
--smoke: development seed 9940, profile D, c 1.0 (runtime check only, not part of any table).
-> results/<job>.json
"""
import os, sys, json, time
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent / "qscan_anchor_theory_20261030" / "scripts"))
import qprofile as qp  # noqa: E402
from env15 import np, lib, semi, sc11, bm, gen, E12OUT  # noqa: E402
from trunc import cut  # noqa: E402

bm.MAXITER = 20000; LAM = 0.5; OUT = HERE / "results"
PROFILES = {"S": [30.0, 200.0, 200.0, 200.0], "D": [250.0, 250.0, 40.0, 40.0],
            "D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0],
            "R1": [100.0, 100.0, 100.0, 100.0], "R2": [60.0, 60.0, 150.0, 150.0], "R3": [150.0, 40.0, 150.0, 150.0],
            "R4": [50.0, 90.0, 250.0, 250.0], "X1": [120.0, 120.0, 70.0, 70.0], "X2": [300.0, 300.0, 60.0, 60.0],
            "X3": [150.0, 150.0, 150.0, 50.0]}
SEEDS = {"S": (9971, 9972, 9973, 9974), "D": (9971, 9972, 9973, 9974), "D2": (9981, 9982, 9983, 9984),
         "M": (9981, 9982, 9983, 9984)}
qp.PROFILES.clear(); qp.PROFILES.update(PROFILES)


def jobs():
    J = [("base", f"T_c{c}_s{s}", None, None, None) for c in (0.65, 0.85, 1.05, 1.25, 1.5) for s in (9701, 9702)]
    for p in PROFILES:
        for s in SEEDS.get(p, (9991, 9992)):
            for c in qp.CS:
                J.append((p, f"c{c}_s{s}", c, s, None))
    for t1 in (2.2, 2.6):                                   # descriptive: D with the analysis window cut (dtrunc data)
        for s in (9971, 9972):
            for c in qp.CS:
                J.append(("D", f"c{c}_s{s}", c, s, t1))
    return J


def anchor_solve(D, cm, pp):
    Q = np.asarray(cm["Q"], float); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * 1.0))
    B, _, scale = lib.dictionary(dict(cm, Q=Q * 1.0), wv["w_in"], wv["w_out"])
    X, info = bm.solve(B, np.eye(D.shape[0]), D, pp["sigma"], scale, pp["K"], LAM, 0.0, pp["live"])
    XL = X[:, pp["live"]]; t = np.asarray(cm["t"], float)
    return dict(t=t.tolist(), a=lib.absorption(cm, t).tolist(), a0=float(lib.absorption(cm, [1.4])[0]),
                e=(XL ** 2).mean(axis=1).tolist(), n=(np.abs(XL) > 1e-8).mean(axis=1).tolist(),
                converged=bool(info["converged"]), iterations=int(info["iterations"]), n_live=int(pp["live"].sum()))


def run(prof, did, c, s, t1, base):
    t0 = time.monotonic()
    if prof == "base":
        cm = semi.curtain_cm(); D = np.load(E12OUT / f"{did}.npz")["Input"].astype(float)
    else:
        ctx = qp.context(base, prof); cm = ctx["cm"]
        g = gen.signal(ctx, c, s); D = g["x"] + qp.MULT * g["level"] * ctx["noise"][0]; del g
    pp = sc11.prepare_sc(D, cm)
    if t1 is not None:
        Dt, ct, pt = cut(D, cm, pp, t1)
        # wavelet from the full section as dtrunc.py; dictionary and solve on the cut window
        Q = np.asarray(cm["Q"], float); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q))
        B, _, scale = lib.dictionary(dict(ct, Q=Q), wv["w_in"], wv["w_out"])
        X, info = bm.solve(B, np.eye(Dt.shape[0]), Dt, pt["sigma"], scale, pt["K"], LAM, 0.0, pt["live"])
        XL = X[:, pt["live"]]; t = np.asarray(ct["t"], float)
        rec = dict(t=t.tolist(), a=lib.absorption(ct, t).tolist(), a0=float(lib.absorption(ct, [1.4])[0]),
                   e=(XL ** 2).mean(axis=1).tolist(), n=(np.abs(XL) > 1e-8).mean(axis=1).tolist(),
                   converged=bool(info["converged"]), iterations=int(info["iterations"]), n_live=int(pt["live"].sum()))
    else:
        rec = anchor_solve(D, cm, pp)
    return dict(profile=prof, id=did, t1=t1, c_anchor=1.0, lam=LAM, seconds=time.monotonic() - t0, **rec)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = gen.prepare(qp.STRUCT, qp.NOISE)
    if "--smoke" in sys.argv:
        r = run("D", "c1.0_s9940", 1.0, 9940, None, base)
        print(r["converged"], r["iterations"], round(r["seconds"]), len(r["t"]), r["n_live"]); sys.exit()
    for prof, did, c, s, t1 in jobs():
        job = f"{prof}_{did}" + (f"_t{t1}" if t1 else ""); path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(prof, did, c, s, t1, base)
        path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, "conv" if rec["converged"] else "NOT-CONVERGED", rec["iterations"], f"{rec['seconds']:.0f}s", flush=True)
