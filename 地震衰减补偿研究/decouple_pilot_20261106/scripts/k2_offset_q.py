"""Pilot K2 (exploratory, no pre-registered rule): can the OFFSET part of the absorption measure Q free of reflectivity colour?
The reflectivity of a CMP is common to its offsets, so ln P_x(f) - ln P_0(f) in a window holds the extra path loss only:
slope = -pi ax / c. The source and the colour cancel. Compared with the vertical estimator (windows of the zero-offset trace,
double centred, same idea as SR), which sees colour trends. c_true in 0.7-1.4, colour trend s (ln amplitude per Hz per s), 20 CMPs
per estimate. Variants: ideal NMO-aligned (the original method's assumption); real NMO stretch on the modelled gather."""
import json
from pathlib import Path
import numpy as np
from core import C, ray, atten, ricker

OUT = Path(__file__).resolve().parents[1] / "results"; OUT.mkdir(exist_ok=True)
N = C.nt; DT = C.dt; T = C.t; OFF = C.offset; CMPS = 20; SEEDS = range(8); WIN = [(a, a + 0.25) for a in np.arange(0.35, 0.96, 0.1)]
BAND = (15.0, 70.0); NF = 1024; f = np.fft.rfftfreq(NF, DT); B = (f >= BAND[0]) & (f <= BAND[1]); MUTE = 1.3
A0 = ray(T, 0)[0]; AX = np.stack([ray(T, x)[1] for x in OFF], axis=1)
vint = np.interp(T, np.r_[0, np.cumsum(2 * C.h / C.vp[:3])], C.vp[:4]); vrms = np.sqrt(np.cumsum(vint ** 2) / np.arange(1, N + 1))
wav = ricker(np.arange(-40, 41) * DT, 45.0)


def reflectivity(rng, s):
    r = (rng.random((N, CMPS)) < 0.3) * rng.standard_normal((N, CMPS))
    if s == 0:
        return r
    L = 100; H = 50; h = 0.5 - 0.5 * np.cos(2 * np.pi * (np.arange(L) + 0.5) / L); ff = np.fft.rfftfreq(2 * L, DT); out = np.zeros_like(r); nm = np.zeros(N)
    for st in range(-H, N, H):
        i0, i1 = max(st, 0), min(st + L, N)
        if i1 - i0 < 2:
            continue
        hh = h[i0 - st:i1 - st]; tc = T[np.clip(st + L // 2, 0, N - 1)]
        out[i0:i1] += np.fft.irfft(np.fft.rfft(r[i0:i1] * hh[:, None], 2 * L, axis=0) * np.exp(s * (tc - 0.6) * (ff - 35))[:, None], 2 * L, axis=0)[:i1 - i0] * hh[:, None]; nm[i0:i1] += hh ** 2
    return out / np.maximum(nm, 1e-6)[:, None]


def conv(x):
    return np.apply_along_axis(lambda v: np.convolve(v, wav, mode="same"), 0, x)


def gather(r, c, x_i, stretch):
    """Attenuated traces of offset x_i, all CMPs. stretch=False: NMO-aligned ideal; True: modelled at the moveout time, then NMO."""
    x = OFF[x_i]
    if not stretch or x == 0:
        return atten((A0 + AX[:, x_i]) / c) @ conv(r)
    tx = np.sqrt(T ** 2 + (x / vrms) ** 2); rp = np.zeros_like(r); pos = tx / DT; i = np.floor(pos).astype(int); fr = pos - i; ok = i + 1 < N
    for j in range(r.shape[1]):
        np.add.at(rp[:, j], i[ok], r[ok, j] * (1 - fr[ok])); np.add.at(rp[:, j], i[ok] + 1, r[ok, j] * fr[ok])
    t0_of = np.interp(T, tx, T); a_pre = (np.interp(t0_of, T, A0) + np.interp(t0_of, T, AX[:, x_i])) / c
    d = atten(a_pre) @ conv(rp); F = np.fft.rfft(d, 4 * N, axis=0); fine = np.fft.irfft(F, 4 * N, axis=0) * 1.0; tf = np.arange(4 * N) * DT / 2
    up = np.fft.irfft(np.fft.rfft(d, axis=0), axis=0, n=d.shape[0])    # keep shape for clarity
    M = 4; Nu = M * N; X = np.fft.rfft(d, axis=0); Xu = np.zeros((Nu // 2 + 1, d.shape[1]), complex); Xu[:X.shape[0]] = X; du = np.fft.irfft(Xu, Nu, axis=0) * M; tu = np.arange(Nu) * DT / M
    return np.stack([np.interp(tx, tu, du[:, j], right=0.0) for j in range(d.shape[1])], axis=1)


def power(y, a, b):
    k = (T >= a) & (T < b); tp = np.hanning(k.sum())[:, None]
    return np.mean(np.abs(np.fft.rfft(y[k] * tp, NF, axis=0)) ** 2, axis=1)


def win_mean(v, a, b):
    k = (T >= a) & (T < b); w = np.hanning(k.sum()) ** 2; return float(w @ v[k] / w.sum())


def smooth(p, m=9):
    return np.convolve(p, np.ones(m) / m, mode="same")


def estimates(r, c, snr, stretch, rng):
    Gs = [gather(r, c, i, stretch) for i in range(len(OFF))]
    sig = np.sqrt(np.mean(Gs[0] ** 2)) / 10 ** (snr / 20); Gn = [g + sig * rng.standard_normal(g.shape) for g in Gs]
    xs, ys = [], []; L0 = []; a0m = []; tp = None
    for a, b in WIN:
        P0 = smooth(power(Gn[0], a, b)); noise = sig ** 2 * np.sum(np.hanning(((T >= a) & (T < b)).sum()) ** 2); Pc0 = np.maximum(P0 - noise, 1e-12 * P0.max())
        L0.append(0.5 * np.log(Pc0[B])); a0m.append(win_mean(A0, a, b)); tc = 0.5 * (a + b)
        for i in range(1, len(OFF)):
            if np.sqrt(tc ** 2 + (OFF[i] / np.interp(tc, T, vrms)) ** 2) / tc > MUTE:
                continue
            Px = smooth(power(Gn[i], a, b)); Pcx = np.maximum(Px - noise, 1e-12 * Px.max())
            sl = np.polyfit(f[B], 0.5 * np.log(Pcx[B] / Pc0[B]), 1)[0]; xs.append(np.pi * win_mean(AX[:, i], a, b)); ys.append(-sl)
    xs, ys = np.array(xs), np.array(ys); c_off = float(np.sum(xs ** 2) / np.sum(xs * ys))
    L0 = np.array(L0); a0m = np.array(a0m); R = L0 - L0.mean(axis=0, keepdims=True); sl = np.array([np.polyfit(f[B], R[k], 1)[0] for k in range(len(WIN))])
    xv = np.pi * (a0m - a0m.mean()); c_ver = float(np.sum(xv ** 2) / np.sum(xv * (-sl)))
    return c_off, c_ver, int(len(xs)), float(np.mean(xs))


rows = []
for stretch in (False, True):
    for snr in (30, 18):
        for s in (0.0, 0.01, -0.01, 0.02, -0.02):
            for c in (0.7, 1.0, 1.4):
                for seed in SEEDS:
                    rng = np.random.default_rng(1000 * seed + 7); r = reflectivity(rng, s)
                    a, b, m, ax = estimates(r, c, snr, stretch, rng)
                    rows.append(dict(stretch=stretch, snr=snr, s=s, c=c, seed=seed, off=np.log(a / c) if a > 0 else None, ver=np.log(b / c) if b > 0 else None, n=m, mean_pi_ax=ax))
            print("stretch", stretch, "snr", snr, "s", s, flush=True)
(OUT / "k2_rows.json").write_text(json.dumps(rows), encoding="utf-8")
print("\nlog(c_est / c_true): mean, sd over c = 0.7, 1.0, 1.4 and 8 seeds (20 CMPs each). off = offset criterion, ver = vertical windows")
for stretch in (False, True):
    for snr in (30, 18):
        print(f"--- {'NMO stretch' if stretch else 'ideal aligned'}, SNR {snr} dB")
        for s in (0.0, 0.01, -0.01, 0.02, -0.02):
            R = [x for x in rows if x["stretch"] == stretch and x["snr"] == snr and x["s"] == s]
            o = np.array([x["off"] for x in R if x["off"] is not None], float); v = np.array([x["ver"] for x in R if x["ver"] is not None], float)
            print(f"s {s:+.3f}: off {np.nanmean(o):+.3f} sd {np.nanstd(o):.3f} ({len(o)}/{len(R)} valid) | ver {np.nanmean(v):+.3f} sd {np.nanstd(v):.3f} ({len(v)}/{len(R)})", flush=True)
print("mean pi*ax per window-offset used:", np.mean([x["mean_pi_ax"] for x in rows]), "number of slopes:", rows[0]["n"], "vertical pi*a0 range:", np.pi * (ray(0.95, 0)[0][0] - ray(0.45, 0)[0][0]))
