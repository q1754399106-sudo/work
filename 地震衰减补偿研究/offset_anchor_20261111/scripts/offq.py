"""Offset spectral ratio on gathers WITHOUT NMO, windows along the moveout curve: absolute cumulative t* (exploratory
synthetic pilot, NOTE_ZH.md). Model and reflectivity generator from decouple_pilot_20261106 (core.py copied unchanged,
k2_offset_q.py); differences: 200 CMPs, events placed at the moveout time by band-limited (Fourier) shifts instead of
linear interpolation, no NMO, adaptive band (power > 3 x noise power, contiguous from 15 Hz), blue reflectivity option,
velocity error in the window positions and in the predicted extra path.

c_off: sum(x^2) / sum(x y), x = pi * extra t* of the offset ray (c = 1), y = -slope of 0.5 ln(P_x / P_0).
c_ver: vertical windows of the zero-offset traces, double centred (as SR; sees colour trends), from k2_offset_q.py.

usage: offq.py run
"""
import sys, json
from pathlib import Path
import numpy as np
from core import C, ray, atten, ricker

OUT = Path(__file__).resolve().parents[1] / "results"; N = C.nt; DT = C.dt; T = C.t; OFF = C.offset; CMPS = 200; SEEDS = range(8)
WIN = [(a, a + 0.25) for a in np.arange(0.35, 0.96, 0.1)]; BAND = (15.0, 70.0); NF = 1024; f = np.fft.rfftfreq(NF, DT); B = (f >= BAND[0]) & (f <= BAND[1]); MUTE = 1.3
vint = np.interp(T, np.r_[0, np.cumsum(2 * C.h / C.vp[:3])], C.vp[:4]); vrms = np.sqrt(np.cumsum(vint ** 2) / np.arange(1, N + 1))
wav = ricker(np.arange(-40, 41) * DT, 45.0); FT = np.fft.rfftfreq(2 * N, DT)


def rays(scale=1.0):
    """a0, ax (n_t, n_off) of the model with all velocities x scale (the processing model when scale != 1)."""
    vp = C.vp; C.vp = vp * scale
    try:
        return ray(T, 0)[0], np.stack([ray(T, x)[1] for x in OFF], axis=1)
    finally:
        C.vp = vp


A0, AX = rays()


def reflectivity(rng, s, beta):
    """White sparse series; optional colour trend s (ln amplitude per Hz per s, as k2) and blue spectrum f^beta."""
    r = (rng.random((N, CMPS)) < 0.3) * rng.standard_normal((N, CMPS))
    if s != 0:
        L = 100; H = 50; h = 0.5 - 0.5 * np.cos(2 * np.pi * (np.arange(L) + 0.5) / L); ff = np.fft.rfftfreq(2 * L, DT); out = np.zeros_like(r); nm = np.zeros(N)
        for st in range(-H, N, H):
            i0, i1 = max(st, 0), min(st + L, N)
            if i1 - i0 < 2:
                continue
            hh = h[i0 - st:i1 - st]; tc = T[np.clip(st + L // 2, 0, N - 1)]
            out[i0:i1] += np.fft.irfft(np.fft.rfft(r[i0:i1] * hh[:, None], 2 * L, axis=0) * np.exp(s * (tc - 0.6) * (ff - 35))[:, None], 2 * L, axis=0)[:i1 - i0] * hh[:, None]; nm[i0:i1] += hh ** 2
        r = out / np.maximum(nm, 1e-6)[:, None]
    if beta:
        r = np.fft.irfft(np.fft.rfft(r, axis=0) * (np.maximum(np.fft.rfftfreq(N, DT), 1.0) ** beta)[:, None], N, axis=0)
    return r


def conv(x):
    return np.apply_along_axis(lambda v: np.convolve(v, wav, mode="same"), 0, x)


def gather(r, c, i):
    """Recorded (no NMO) traces of offset i, all CMPs: event of t0 at t_x = sqrt(t0^2 + (x / vrms)^2), band-limited shift."""
    if OFF[i] == 0:
        return atten(A0 / c) @ conv(r)
    tx = np.sqrt(T ** 2 + (OFF[i] / vrms) ** 2); ok = tx < T[-1] - 0.1
    rp = np.fft.irfft(np.exp(-2j * np.pi * FT[:, None] * tx[None, ok]) @ r[ok], 2 * N, axis=0)[:N]
    t0 = np.interp(T, tx, T); a = (np.interp(t0, T, A0) + np.interp(t0, T, AX[:, i])) / c
    return atten(a) @ conv(rp)


def power(y, a, b):
    k = (T >= a) & (T < b); h = np.hanning(k.sum())
    return np.convolve(np.mean(np.abs(np.fft.rfft(y[k] * h[:, None], NF, axis=0)) ** 2, axis=1), np.ones(9) / 9, mode="same"), np.sum(h ** 2)


def wmean(v, a, b):
    k = (T >= a) & (T < b); w = np.hanning(k.sum()) ** 2; return float(w @ v[k] / w.sum())


def estimates(r, c, snr, vscale, rng):
    G = [gather(r, c, i) for i in range(len(OFF))]; sig = np.sqrt(np.mean(G[0] ** 2)) / 10 ** (snr / 20) if snr else 0.0
    G = [g + sig * rng.standard_normal(g.shape) for g in G]; a0p, axp = (A0, AX) if vscale == 1 else rays(vscale); vp = vrms * vscale
    xs, ys, L0, a0m = [], [], [], []
    for a, b in WIN:
        R0, e = power(G[0], a, b); n0 = sig ** 2 * e; P0 = np.maximum(R0 - n0, 1e-12 * R0.max()); L0.append(0.5 * np.log(P0[B])); a0m.append(wmean(a0p, a, b)); tc = 0.5 * (a + b)
        for i in range(1, len(OFF)):
            ta, tm, tb = np.sqrt(np.array([a, tc, b]) ** 2 + (OFF[i] / np.interp([a, tc, b], T, vp)) ** 2)
            if tb > T[-1] - 0.1 or tm / tc > MUTE:
                continue
            Rx, ex = power(G[i], ta, tb); nx = sig ** 2 * ex; Px = np.maximum(Rx - nx, 1e-12 * Rx.max())
            ok = B & (R0 > 3 * n0) & (Rx > 3 * nx); ok &= f < (f[B][np.argmin(ok[B])] if not ok[B].all() else 1e9)
            if ok.sum() < 12:
                continue
            xs.append(np.pi * wmean(axp[:, i], a, b)); ys.append(-np.polyfit(f[ok], 0.5 * np.log(Px[ok] / P0[ok]), 1)[0])
    xs, ys = np.array(xs), np.array(ys); c_off = float(np.sum(xs ** 2) / np.sum(xs * ys)) if len(xs) else float("nan")
    L0 = np.array(L0); R = L0 - L0.mean(axis=0, keepdims=True); sl = np.array([np.polyfit(f[B], R[k], 1)[0] for k in range(len(WIN))])
    xv = np.pi * (np.array(a0m) - np.mean(a0m)); c_ver = float(np.sum(xv ** 2) / np.sum(xv * (-sl)))
    return c_off, c_ver, len(xs)


CASES = [dict(tag="noise-free white", snr=0, s=0.0, beta=0.0, v=1.0)] + \
        [dict(tag="white", snr=q, s=0.0, beta=0.0, v=1.0) for q in (30, 18, 10)] + \
        [dict(tag="blue", snr=q, s=0.0, beta=0.7, v=1.0) for q in (30, 18, 10)] + \
        [dict(tag="blue, trend %+.2f" % s, snr=18, s=s, beta=0.7, v=1.0) for s in (0.02, -0.02)] + \
        [dict(tag="blue, velocity %+d%%" % round(100 * (v - 1)), snr=18, s=0.0, beta=0.7, v=v) for v in (1.02, 0.98, 1.05, 0.95)]

if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); rows = []
    for case in CASES:
        for c in (0.7, 1.0, 1.4):
            for seed in SEEDS:
                rng = np.random.default_rng(1000 * seed + 7); r = reflectivity(rng, case["s"], case["beta"])
                co, cv, n = estimates(r, c, case["snr"], case["v"], rng); rows.append(dict(case, c=c, seed=seed, c_off=co, c_ver=cv, n=n))
            g = [q for q in rows if q["tag"] == case["tag"] and q["snr"] == case["snr"] and q["c"] == c]
            lo = np.array([np.log(q["c_off"] / c) if q["c_off"] > 0 else np.nan for q in g]); lv = np.array([np.log(q["c_ver"] / c) if q["c_ver"] > 0 else np.nan for q in g])
            print("%-22s snr %2d c %.1f | offset %+.3f sd %.3f (bad %d, n %d) | vertical %+.3f sd %.3f (bad %d)" % (case["tag"], case["snr"], c, np.nanmean(lo), np.nanstd(lo), np.isnan(lo).sum(), g[0]["n"], np.nanmean(lv), np.nanstd(lv), np.isnan(lv).sum()), flush=True)
    (OUT / "offq_rows.json").write_text(json.dumps(rows), encoding="utf-8")
