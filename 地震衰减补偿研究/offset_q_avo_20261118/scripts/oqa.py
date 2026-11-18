"""Offset spectral ratio on gathers WITHOUT NMO (windows along the moveout curve) with Shuey AVO in the reflectivity:
how well does it measure c = Q_true / Q_model (one global factor) when N CMPs are pooled? (PLAN_GATE_A1_ZH.md)
Estimator from offset_anchor_20261111/scripts/offq.py (c_off), changed: model Q per layer, Ricker 40/60 Hz, band 12-90 Hz
(adaptive: power > 3 x noise, contiguous from the low edge), Shuey AVO r(t0, x) = A + B sin^2(theta), per-trace noise
(rms of each offset trace, as step 2b), estimates on groups of N CMPs out of 200.
usage: oqa.py run <model> [dev] | smoke"""
import sys, json, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import C, ray, atten, ricker  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; N = C.nt; DT = C.dt; T = C.t; OFF = C.offset; NX = len(OFF); CMPS = 200
WIN = [(a, a + 0.25) for a in np.arange(0.35, 0.96, 0.1)]; BAND = (12.0, 90.0); NF = 1024; f = np.fft.rfftfreq(NF, DT)
B = (f >= BAND[0]) & (f <= BAND[1]); MUTE = 1.3; FT = np.fft.rfftfreq(2 * N, DT)
vint = np.interp(T, np.r_[0, np.cumsum(2 * C.h / C.vp[:3])], C.vp[:4]); vrms = np.sqrt(np.cumsum(vint ** 2) / np.arange(1, N + 1))
MODELS = {"M1": [40, 50, 70, 100], "M16": [45, 60, 35, 80], "M17": [90, 35, 55, 45], "M18": [32, 50, 90, 60], "M19": [65, 85, 30, 40]}
CS = (0.7, 1.0, 1.4); FPS = (40.0, 60.0); SEEDS = range(8); GROUPS = (10, 50, 200)
CASES = [("avo", 24, 1.0), ("avo", 30, 1.0), ("avo", 18, 1.0), ("flat", 24, 1.0), ("avo", 24, 1.02), ("avo", 24, 0.98)]


def rays(Q, scale=1.0):
    vp = C.vp; C.vp = vp * scale
    try:
        a0 = ray(T, 0, Q)[0]; out = [ray(T, x, Q) for x in OFF]
        return a0, np.stack([o[1] for o in out], 1), np.stack([o[2] for o in out], 1)
    finally:
        C.vp = vp


def reflectivity(rng, avo):
    """A: sparse (30 %) standard normal; B = -0.8 A + 0.6 e (e ~ N(0,1) on the same spikes), as the step-2b dense scene ratio."""
    on = rng.random((N, CMPS)) < 0.3; A = on * rng.standard_normal((N, CMPS)); Bg = (-0.8 * A + 0.6 * on * rng.standard_normal((N, CMPS))) if avo else 0 * A
    return A, Bg


def gathers(A, Bg, Q, c, fp, P):
    """Noise-free recorded traces per offset (N, CMPS). True path t* = model path / c. P: ray parameters (N, NX) at t0."""
    a0, ax, _ = rays(Q); wav = ricker(np.arange(-40, 41) * DT, fp); out = []
    conv = lambda x: np.apply_along_axis(lambda v: np.convolve(v, wav, mode="same"), 0, x)
    for i in range(NX):
        s2 = (P[:, i] * vint) ** 2; r = A + Bg * s2[:, None]
        if OFF[i] == 0:
            out.append(atten(a0 / c) @ conv(r)); continue
        tx = np.sqrt(T ** 2 + (OFF[i] / vrms) ** 2); ok = tx < T[-1] - 0.1
        rp = np.fft.irfft(np.exp(-2j * np.pi * FT[:, None] * tx[None, ok]) @ r[ok], 2 * N, axis=0)[:N]
        t0 = np.interp(T, tx, T); a = (np.interp(t0, T, a0) + np.interp(t0, T, ax[:, i])) / c
        out.append(atten(a) @ conv(rp))
    return out


def power(y, a, b):
    k = (T >= a) & (T < b); h = np.hanning(k.sum())
    return np.convolve(np.mean(np.abs(np.fft.rfft(y[k] * h[:, None], NF, axis=0)) ** 2, axis=1), np.ones(9) / 9, mode="same"), np.sum(h ** 2)


def wmean(v, a, b):
    k = (T >= a) & (T < b); w = np.hanning(k.sum()) ** 2; return float(w @ v[k] / w.sum())


def c_off(G, sig, a0p, axp, vp):
    xs, ys = [], []
    for a, b in WIN:
        R0, e = power(G[0], a, b); n0 = sig[0] ** 2 * e; P0 = np.maximum(R0 - n0, 1e-12 * R0.max()); tc = 0.5 * (a + b)
        for i in range(1, NX):
            ta, tm, tb = np.sqrt(np.array([a, tc, b]) ** 2 + (OFF[i] / np.interp([a, tc, b], T, vp)) ** 2)
            if tb > T[-1] - 0.1 or tm / tc > MUTE:
                continue
            Rx, ex = power(G[i], ta, tb); nx = sig[i] ** 2 * ex; Px = np.maximum(Rx - nx, 1e-12 * Rx.max())
            ok = B & (R0 > 3 * n0) & (Rx > 3 * nx); ok &= f < (f[B][np.argmin(ok[B])] if not ok[B].all() else 1e9)
            if ok.sum() < 12:
                continue
            xs.append(np.pi * wmean(axp[:, i], a, b)); ys.append(-np.polyfit(f[ok], 0.5 * np.log(Px[ok] / P0[ok]), 1)[0])
    xs, ys = np.array(xs), np.array(ys)
    return (float(np.sum(xs ** 2) / np.sum(xs * ys)) if len(xs) and np.sum(xs * ys) > 0 else float("nan")), len(xs)


def run(model, dev=False):
    Q = np.array(MODELS[model], float); a0m, axm, P = rays(Q); rows = []; t0 = time.monotonic(); mi = list(MODELS).index(model)
    proc = {s: rays(Q, s)[:2] for s in (1.0, 1.02, 0.98)}
    seeds = range(2) if dev else SEEDS; cs = (1.0,) if dev else CS
    for fp in FPS:
        for c in cs:
            for seed in seeds:
                for avo in ("avo", "flat"):
                    rng = np.random.default_rng(31_000_000 + mi * 100_000 + int(fp) * 1000 + int(c * 10) * 10 + seed)
                    A, Bg = reflectivity(rng, avo == "avo"); G0 = gathers(A, Bg, Q, c, fp, P); rms = np.array([np.sqrt(np.mean(g ** 2)) for g in G0])
                    for kind, snr, vs in CASES:
                        if kind != avo:
                            continue
                        nrng = np.random.default_rng(rng.integers(1 << 62) + snr); sig = rms / 10 ** (snr / 20)
                        G = [g + s * nrng.standard_normal(g.shape) for g, s in zip(G0, sig)]; a0p, axp = proc[vs]
                        for n in GROUPS:
                            for k in range(CMPS // n):
                                cols = slice(k * n, (k + 1) * n); ce, nw = c_off([g[:, cols] for g in G], sig, a0p, axp, vrms * vs)
                                rows.append(dict(model=model, fp=fp, c=c, seed=seed, kind=kind, snr=snr, v=vs, n=n, group=k, c_hat=ce, nwin=nw))
                print(model, fp, c, seed, "%.0fs" % (time.monotonic() - t0), flush=True)
    (ROOT / "results" / ("dev_%s.json" % model if dev else "a1_%s.json" % model)).write_text(json.dumps(rows))


def smoke():
    """Noise-free, AVO off, c = 1, all 200 CMPs, used model M1: the estimator runs and returns a finite value (no score printed)."""
    Q = np.array(MODELS["M1"], float); a0m, axm, P = rays(Q); rng = np.random.default_rng(5); A, Bg = reflectivity(rng, False)
    G = gathers(A, Bg, Q, 1.0, 40.0, P); ce, nw = c_off(G, np.zeros(NX), a0m, axm, vrms); print("finite", np.isfinite(ce), "windows", nw)


if __name__ == "__main__":
    if sys.argv[1] == "smoke":
        smoke()
    else:
        run(sys.argv[2], dev=len(sys.argv) > 3)
