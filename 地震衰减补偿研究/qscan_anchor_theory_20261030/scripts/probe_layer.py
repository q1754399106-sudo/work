"""Exploratory probe (THEORY_ZH.md section 5i): two-layer Q scan, deep-layer factor only.

usage: probe_layer.py    prints lambda, c_true, v2, objective, converged for development seed 9940.
Experiment-14 recipe (structure line 10227, noise x 0.5), wavelet at the correct anchor (c_a = 1, lib.wavelets),
scanned model Q / [u, u, v2, v2] (shallow factor fixed at the truth), lasso objective of lib.qscan.
"""
from env15 import np, lib, semi, sc11, bm, gen
bm.MAXITER = 20000
ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
for lam, cs in ((0.5, (0.7, 1.0, 1.4)), (2 ** -4, (1.0,))):
    for C in cs:
        g = gen.signal(ctx, C, 9940); D = g["x"] + 0.5 * g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
        wv = lib.wavelets(D, semi.TW, cm); n = D.shape[0]; sigma, K, live = pp["sigma"], pp["K"], pp["live"]
        s2 = np.where(live, sigma ** 2, np.inf); prev = None
        for v2 in np.arange(1.6, -0.61, -0.2):
            with np.errstate(divide="ignore"):
                qv = Q / np.array([1 / C, 1 / C, v2, v2])
            B, _, scale = lib.dictionary(dict(cm, Q=qv), wv["w_in"], wv["w_out"])
            pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
            X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, lam, 0.0, live, start=prev); prev = X; res = D - B @ X
            J = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + lam * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
            print(lam, C, round(float(v2), 2), round(J, 1), info["converged"], flush=True)
