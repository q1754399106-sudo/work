"""Data for paper figure 5 (diag_bias.py logic, saved): noise-free M1, AVO off, c = 1, 200 CMPs, 4 seeds. ln(c_hat/c) vs wavelet
peak frequency for the original estimator and the equal-length, unsmoothed one; per-window y/x vs offset at 40 Hz."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import oqa as o

T = o.T; f = o.f


def power(y, a, b, smooth, L=None):
    if L is not None:
        m = 0.5 * (a + b); a, b = m - L / 2, m + L / 2
    k = (T >= a) & (T < b); h = np.hanning(k.sum()); P = np.mean(np.abs(np.fft.rfft(y[k] * h[:, None], o.NF, axis=0)) ** 2, axis=1)
    return np.convolve(P, np.ones(smooth) / smooth, mode="same") if smooth > 1 else P


def c_off(G, a0p, axp, smooth, eqlen):
    B = (f >= 12.0) & (f <= 90.0); xs, ys, per = [], [], []
    for a, b in o.WIN:
        P0 = power(G[0], a, b, smooth); tc = 0.5 * (a + b)
        for i in range(1, o.NX):
            ta, tm, tb = np.sqrt(np.array([a, tc, b]) ** 2 + (o.OFF[i] / np.interp([a, tc, b], T, o.vrms)) ** 2)
            if tb > T[-1] - 0.1 or tm / tc > o.MUTE:
                continue
            Px = power(G[i], ta, tb, smooth, L=(b - a) if eqlen else None); ok = B & (P0 > 1e-6 * P0.max()) & (Px > 1e-6 * Px.max())
            if ok.sum() < 12:
                continue
            x = np.pi * o.wmean(axp[:, i], a, b); y = -np.polyfit(f[ok], 0.5 * np.log(Px[ok] / P0[ok]), 1)[0]; xs.append(x); ys.append(y); per.append((float(tc), float(o.OFF[i]), float(y / x)))
    xs, ys = np.array(xs), np.array(ys); return float(np.sum(xs ** 2) / np.sum(xs * ys)), per


Q = np.array(o.MODELS["M1"], float); a0m, axm, P = o.rays(Q); out = dict(bias=[], per_offset=[])
for fp in (30.0, 35.0, 40.0, 50.0, 60.0, 80.0):
    for seed in range(4):
        rng = np.random.default_rng(777 + seed); A, Bg = o.reflectivity(rng, False); G = o.gathers(A, Bg, Q, 1.0, fp, P)
        for name, sm, eq in (("original", 9, False), ("equal_length", 1, True)):
            ce, per = c_off(G, a0m, axm, sm, eq); out["bias"].append(dict(fp=fp, seed=seed, est=name, ln_err=float(np.log(ce))))
            if fp == 40.0:
                out["per_offset"] += [dict(seed=seed, est=name, tc=p[0], offset=p[1], y_over_x=p[2]) for p in per]
    print(fp, flush=True)
(o.ROOT / "results" / "fig_window.json").write_text(json.dumps(out))
