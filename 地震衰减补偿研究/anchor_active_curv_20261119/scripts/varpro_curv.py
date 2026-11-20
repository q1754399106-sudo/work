"""PLAN_VARPRO_ZH.md: anchor-point lasso solve (same setting as active_curv.py) plus
  VP  variable-projection Gauss-Newton curvature  Phi = (P_perp Q)^T (P_perp Q),  q_k = X_kj d b_k / d A_k  (active set per trace)
  NF  non-zero density / atom spectral centroid
  G   diagonal, no projection: sum_j X_kj^2 |d b_k / d A_k|^2   (descriptive)
  E, N reproduction of the archived energy envelope and non-zero density.
usage: varpro_curv.py [--smoke]   worker; jobs claimed with results_vp/<job>.claim, resumable.
Data and caches are read from the main checkout (MAIN); results are written next to this script's experiment folder.
"""
import os, sys, json, time
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
MAIN = Path(r"D:\ChatGPT\研究\地震衰减补偿研究")
sys.path.insert(0, str(MAIN / "qscan_anchor_theory_20261030" / "scripts"))
import qprofile as qp  # noqa: E402
from env15 import np, lib, semi, sc11, bm, gen, E12OUT  # noqa: E402

bm.MAXITER = 20000; LAM = 0.5; OUT = HERE / "results_vp"; H = 1e-4
PROFILES = {"S": [30.0, 200.0, 200.0, 200.0], "D": [250.0, 250.0, 40.0, 40.0],
            "D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0],
            "R1": [100.0, 100.0, 100.0, 100.0], "R2": [60.0, 60.0, 150.0, 150.0], "R3": [150.0, 40.0, 150.0, 150.0],
            "R4": [50.0, 90.0, 250.0, 250.0], "X1": [120.0, 120.0, 70.0, 70.0], "X2": [300.0, 300.0, 60.0, 60.0],
            "X3": [150.0, 150.0, 150.0, 50.0]}
SEEDS = {"S": (9971, 9972, 9973, 9974), "D": (9971, 9972, 9973, 9974), "D2": (9981, 9982, 9983, 9984),
         "M": (9981, 9982, 9983, 9984)}
qp.PROFILES.clear(); qp.PROFILES.update(PROFILES)


def jobs():
    J = [("base", f"T_c{c}_s{s}", None, None) for c in (0.65, 0.85, 1.05, 1.25, 1.5) for s in (9701, 9702)]
    for p in PROFILES:
        for s in SEEDS.get(p, (9991, 9992)):
            for c in qp.CS:
                J.append((p, f"c{c}_s{s}", c, s))
    return J


def dict_a(cm, w_in, a):
    """lib.dictionary (normalised B only) with the per-atom two-way t* given explicitly."""
    t = bm.vec(cm["t"]); n = len(t); N = int(cm["nfft"]); f = np.arange(N // 2 + 1) / (N * float(cm["dt"]))
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); b[0] = 0
    G = np.empty((n, n)); ph = np.exp(-2j * np.pi * np.arange(N // 2 + 1)[:, None] * np.arange(n)[None, :] / N)
    for k in range(n):
        Hk = np.exp(b * a[k]); Hk[-1] = Hk[-1].real
        G[:, k] = np.fft.irfft(Hk * ph[:, k], n=N)[:n]
    B = G @ lib.toeplitz_wavelet(w_in, n)
    return B / np.linalg.norm(B, axis=0)


def anchor(D, cm, pp):
    Q = np.asarray(cm["Q"], float); cmq = dict(cm, Q=Q * 1.0); wv = lib.wavelets(D, semi.TW, cmq)
    B, _, scale = lib.dictionary(cmq, wv["w_in"], wv["w_out"])
    X, info = bm.solve(B, np.eye(D.shape[0]), D, pp["sigma"], scale, pp["K"], LAM, 0.0, pp["live"])
    t = np.asarray(cm["t"], float); a = lib.absorption(cmq, t); a0 = float(lib.absorption(cmq, [1.4])[0])
    chk = float(np.max(np.abs(dict_a(cmq, wv["w_in"], a) - B)))
    dB = (dict_a(cmq, wv["w_in"], a + H) - dict_a(cmq, wv["w_in"], a - H)) / (2 * H)
    live = np.flatnonzero(pp["live"]); XL = X[:, live]
    num = den = trphi = tot = 0.0; phiG = np.zeros(len(t)); nact = []
    dn2 = np.sum(dB ** 2, axis=0)
    for j in live:
        idx = np.flatnonzero(X[:, j] != 0); nact.append(len(idx))
        if len(idx) == 0:
            continue
        DA = B[:, idx]; Qm = dB[:, idx] * X[idx, j]
        R = Qm - DA @ np.linalg.lstsq(DA, Qm, rcond=None)[0]
        va = R @ a[idx]; v1 = R.sum(axis=1)
        num += float(va @ v1); den += float(va @ va); trphi += float(np.sum(R * R)); tot += float(v1 @ v1)
        phiG[idx] += X[idx, j] ** 2 * dn2[idx]
    S = np.abs(np.fft.rfft(B, 1024, axis=0)) ** 2; fr = np.fft.rfftfreq(1024, float(cm["dt"]))
    fc = (fr @ S) / S.sum(axis=0)
    return dict(t=t.tolist(), a=a.tolist(), a0=a0, e=(XL ** 2).mean(axis=1).tolist(), n=(np.abs(XL) > 1e-8).mean(axis=1).tolist(),
                num_vp=num, den_vp=den, alpha_vp=a0 * num / den, ratio_1phi1_tr=tot / trphi, phiG=phiG.tolist(), fc=fc.tolist(),
                n_active_median=float(np.median(nact)), dict_check=chk,
                converged=bool(info["converged"]), iterations=int(info["iterations"]), n_live=int(pp["live"].sum()))


def run(prof, did, c, s, base):
    t0 = time.monotonic()
    if prof == "base":
        cm = semi.curtain_cm(); D = np.load(E12OUT / f"{did}.npz")["Input"].astype(float)
    else:
        ctx = qp.context(base, prof); cm = ctx["cm"]
        g = gen.signal(ctx, c, s); D = g["x"] + qp.MULT * g["level"] * ctx["noise"][0]; del g
    rec = anchor(D, cm, sc11.prepare_sc(D, cm))
    return dict(profile=prof, id=did, c_anchor=1.0, lam=LAM, h=H, seconds=time.monotonic() - t0, **rec)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = gen.prepare(qp.STRUCT, qp.NOISE)
    if "--smoke" in sys.argv:
        r = run("D", "c1.0_s9940", 1.0, 9940, base)
        print(r["converged"], r["iterations"], round(r["seconds"]), r["n_live"], "alpha_vp %.3f" % r["alpha_vp"],
              "ratio %.2f" % r["ratio_1phi1_tr"], "dict_check %.2e" % r["dict_check"], "n_act %.0f" % r["n_active_median"]); sys.exit()
    for prof, did, c, s in jobs():
        job = f"{prof}_{did}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(prof, did, c, s, base)
        path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, "conv" if rec["converged"] else "NOT-CONVERGED", rec["iterations"], f"{rec['seconds']:.0f}s",
              "vp %.3f" % rec["alpha_vp"], flush=True)
