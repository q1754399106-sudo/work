"""Development only (used model M1, noise-free, AVO off, c = 1, 200 CMPs, 4 seeds): which estimator setting carries
the wavelet-dependent bias of c_off? Variants change one thing at a time."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import oqa as o

T = o.T; f = o.f


def power(y, a, b, smooth=9, L=None):
    if L is not None:  # fixed window length centred in [a, b]
        m = 0.5 * (a + b); a, b = m - L / 2, m + L / 2
    k = (T >= a) & (T < b); h = np.hanning(k.sum())
    P = np.mean(np.abs(np.fft.rfft(y[k] * h[:, None], o.NF, axis=0)) ** 2, axis=1)
    return (np.convolve(P, np.ones(smooth) / smooth, mode="same") if smooth > 1 else P)


def c_off(G, a0p, axp, smooth=9, eqlen=False, band=(12.0, 90.0), xmode="hann2", perwin=False):
    B = (f >= band[0]) & (f <= band[1]); xs, ys, per = [], [], []
    for a, b in o.WIN:
        P0 = power(G[0], a, b, smooth); tc = 0.5 * (a + b)
        for i in range(1, o.NX):
            ta, tm, tb = np.sqrt(np.array([a, tc, b]) ** 2 + (o.OFF[i] / np.interp([a, tc, b], T, o.vrms)) ** 2)
            if tb > T[-1] - 0.1 or tm / tc > o.MUTE:
                continue
            Px = power(G[i], ta, tb, smooth, L=(b - a) if eqlen else None)
            ok = B & (P0 > 1e-6 * P0.max()) & (Px > 1e-6 * Px.max())
            if ok.sum() < 12:
                continue
            if xmode == "hann2":
                x = np.pi * o.wmean(axp[:, i], a, b)
            else:  # centre value
                x = np.pi * np.interp(tc, T, axp[:, i])
            y = -np.polyfit(f[ok], 0.5 * np.log(Px[ok] / P0[ok]), 1)[0]; xs.append(x); ys.append(y); per.append((tc, i, y / x))
    xs, ys = np.array(xs), np.array(ys)
    return float(np.sum(xs ** 2) / np.sum(xs * ys)), per


Q = np.array(o.MODELS["M1"], float); a0m, axm, P = o.rays(Q)
VARIANTS = {"base": {}, "no_smooth": dict(smooth=1), "eq_len": dict(eqlen=True), "band20-80": dict(band=(20.0, 80.0)),
            "band12-60": dict(band=(12.0, 60.0)), "x_centre": dict(xmode="centre"), "eq_len+no_smooth": dict(eqlen=True, smooth=1)}
for fp in (30.0, 40.0, 50.0, 60.0, 80.0):
    res = {k: [] for k in VARIANTS}; per_all = []
    for seed in range(4):
        rng = np.random.default_rng(777 + seed); A, Bg = o.reflectivity(rng, False); G = o.gathers(A, Bg, Q, 1.0, fp, P)
        for k, kw in VARIANTS.items():
            ce, per = c_off(G, a0m, axm, **kw); res[k].append(np.log(ce))
            if k == "base":
                per_all += per
    print("fp %4.0f  " % fp + "  ".join("%s %+.3f" % (k, np.mean(v)) for k, v in res.items()), flush=True)
    per_all = np.array(per_all)
    for i in (2, 5, 8, 10):
        m = per_all[:, 1] == i
        if m.any():
            print("     base: offset %4.0f m  median y/x %.3f (n %d)" % (o.OFF[i], np.median(per_all[m, 2]), m.sum()))
    for tc in np.unique(per_all[:, 0]):
        m = per_all[:, 0] == tc
        print("     base: window centre %.3f s  median y/x %.3f" % (tc, np.median(per_all[m, 2])))
