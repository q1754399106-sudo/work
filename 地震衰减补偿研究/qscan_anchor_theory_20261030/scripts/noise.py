"""Noise test (PLAN_NOISE_ZH.md): shift invariance of the Q selection against the noise level.

usage: noise.py    worker; (dataset, scan) jobs claimed with results3/<job>.claim, resumable.
Datasets: experiment-14 recipe (gen.signal on structure line 10227), c 0.7 / 1.0 / 1.4, seeds 9931 / 9932,
noise x 0.25 / 0.5 / 1 / 2. Scans: QScan-W at lambda 1/16; fixed-anchor QScan at lambda 1/2 with the wavelet window
1.2-1.6 s or 2.4-2.8 s and c_a 0.7 or 1.4.
-> results3/<dataset>__<scan>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from anchor2 import wavelets_at

STRUCT, NOISE = 10227, [10149, 10273, 10355]
CS = [0.7, 1.0, 1.4]; SEEDS = [9931, 9932]; MULTS = [0.25, 0.5, 1.0, 2.0]
WINDOWS = {"w1.4": (1.2, 1.6, 1.4), "w2.6": (2.4, 2.8, 2.6)}; ANCHORS = [0.7, 1.4]
LAM_W, LAM_A, MAXITER = 2 ** -4, 0.5, 20000
OUT = ROOT / "results3"; bm.MAXITER = MAXITER


def scans():
    return ["qscanw"] + [f"{w}_a{a}" for w in WINDOWS for a in ANCHORS]


def jobs():
    return [(f"m{m}_c{c}_s{s}", m, c, s, sc) for m in MULTS for c in CS for s in SEEDS for sc in scans()]


def scan(D, cm, Q, pp, lam, wave):
    """Lasso objective of lib.qscan over QGRID; wave(cc) returns the wavelets for scanned model cc."""
    sigma, K, live = pp["sigma"], pp["K"], pp["live"]; n = D.shape[0]
    J = np.zeros(len(lib.QGRID)); C = np.zeros(J.shape, bool); IT = np.zeros(J.shape, int)
    s2 = np.where(live, sigma ** 2, np.inf); prev = None
    for j, c in enumerate(lib.QGRID):
        cc = dict(cm); cc["Q"] = Q * c; wv = wave(cc)
        B, _, scale = lib.dictionary(cc, wv["w_in"], wv["w_out"])
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live, start=prev); prev = X
        res = D - B @ X
        J[j] = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + lam * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
        C[j], IT[j] = bool(info["converged"]), int(info["iterations"])
    return float(lib.QGRID[int(np.argmin(J))]), dict(J=J.tolist(), converged=C.tolist(), iterations=IT.tolist())


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(STRUCT, NOISE); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    for did, m, c, s, sc in jobs():
        job = f"{did}__{sc}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic(); g = gen.signal(ctx, c, s); D = g["x"] + m * g["level"] * ctx["noise"][0]
        pp = sc11.prepare_sc(D, cm)
        if sc == "qscanw":
            lam, ca, win = LAM_W, None, None
            sel, diag = scan(D, cm, Q, pp, lam, lambda cc: lib.wavelets(D, semi.TW, cc))
        else:
            w, a = sc.split("_a"); lam, ca, win = LAM_A, float(a), WINDOWS[w]
            ca_cm = dict(cm); ca_cm["Q"] = Q * ca; wv = wavelets_at(D, semi.TW, ca_cm, *win)
            sel, diag = scan(D, cm, Q, pp, lam, lambda cc: wv)
        rec = dict(id=did, scan=sc, noise_mult=m, c_true=c, seed=s, level=g["level"], lam=lam, c_anchor=ca,
                   window=list(win) if win else None, c_sel=sel, diag=diag, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, sel, f"{rec['seconds']:.0f}s", flush=True)
