"""PLAN_DIAG_RESID_ZH.md: where does the residual positive bias of the frozen estimator (far-offset equal-length window,
no smoothing) come from? Noise-free M1, c = 1, AVO off, seeds 777-780, 200 CMPs, wavelets 30/40/60 Hz. Variants change one
thing per (window, offset) synthesis: V0 as is, V1 no moveout compression, V3 constant attenuation in the window, V2 one
event per window, V2n one event without dispersion phase, V0b (descriptive) fit band limited to -40 dB of each peak.
-> results/diag_resid.json
"""
import sys, json, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import oqa as o  # noqa: E402
from core import C, atten, ricker  # noqa: E402

o.EQLEN = True; o.SMOOTH = 1                       # frozen estimator (DIAG_A1_ZH.md section 3)
T, N, DT, FT, f, OFF = o.T, o.N, o.DT, o.FT, o.f, o.OFF
Q = np.array(o.MODELS["M1"], float)
SEEDS = range(4); FPS = (30.0, 40.0, 60.0); VARS = ("V0", "V0b", "V1", "V3", "V2", "V2n")


def atten_amp(a):
    """As core.atten, amplitude only (no dispersion phase)."""
    M = C.nfft; fr = np.arange(M // 2 + 1) / (M * C.dt); H = np.exp(np.outer(-np.pi * fr, a)); H[-1] = H[-1].real
    return np.fft.irfft(H * np.fft.rfft(np.eye(C.nt), M, axis=0), M, axis=0)[:C.nt]


def moveout(r, tx):
    ok = tx < T[-1] - 0.1
    return np.fft.irfft(np.exp(-2j * np.pi * FT[:, None] * tx[None, ok]) @ r[ok], 2 * N, axis=0)[:N]


def slope(y0, yx, a, b, ta, tb, cut40=False):
    R0, _ = o.power(y0, a, b); Rx, _ = o.power(yx, ta, tb)
    ok = o.B & (R0 > 0) & (Rx > 0)
    if cut40:
        ok &= (R0 > 1e-4 * R0.max()) & (Rx > 1e-4 * Rx.max())
    ok &= f < (f[o.B][np.argmin(ok[o.B])] if not ok[o.B].all() else 1e9)
    if ok.sum() < 12:
        return None
    return -np.polyfit(f[ok], 0.5 * np.log(Rx[ok] / R0[ok]), 1)[0]


def one(fp):
    a0m, axm, P = o.rays(Q); wav = ricker(np.arange(-40, 41) * DT, fp)
    conv = lambda x: np.apply_along_axis(lambda v: np.convolve(v, wav, mode="same"), 0, x)
    out = {v: [] for v in VARS}
    for seed in SEEDS:
        rng = np.random.default_rng(777 + seed); A, Bg = o.reflectivity(rng, False)
        G = o.gathers(A, Bg, Q, 1.0, fp, P); cA = conv(A)
        rng1 = np.random.default_rng(90_777 + seed); amp1 = rng1.standard_normal(o.CMPS)
        acc = {v: ([], []) for v in VARS}
        for a, b in o.WIN:
            tc = 0.5 * (a + b); L = b - a; k1 = int(round(tc / DT)); t1 = T[k1]
            r1 = np.zeros((N, o.CMPS)); r1[k1] = amp1; cr1 = conv(r1)
            a0c = float(np.interp(tc, T, a0m)); a01 = float(a0m[k1])
            near3 = atten(np.full(N, a0c)) @ cA; near2 = atten(np.full(N, a01)) @ cr1; near2n = atten_amp(np.full(N, a01)) @ cr1
            for i in range(1, o.NX):
                ta, tm, tb = np.sqrt(np.array([a, tc, b]) ** 2 + (OFF[i] / np.interp([a, tc, b], T, o.vrms)) ** 2)
                if tb > T[-1] - 0.1 or tm / tc > o.MUTE:
                    continue
                ta, tb = tm - L / 2, tm + L / 2; tx = np.sqrt(T ** 2 + (OFF[i] / o.vrms) ** 2)
                x0 = np.pi * o.wmean(axm[:, i], a, b); axc = float(np.interp(tc, T, axm[:, i])); ax1 = float(axm[k1, i])
                # V0 / V0b
                for v, cut in (("V0", False), ("V0b", True)):
                    y = slope(G[0], G[i], a, b, ta, tb, cut)
                    if y is not None:
                        acc[v][0].append(x0); acc[v][1].append(y)
                # V1: whole window shifted by the centre moveout, relative spacing kept, same total loss per event
                d = tm - tc; tt = T - d
                arec = np.interp(tt, T, a0m) + np.interp(tt, T, axm[:, i])
                far1 = atten(arec) @ conv(moveout(A, T + d))
                y = slope(G[0], far1, a, b, ta, tb)
                if y is not None:
                    acc["V1"][0].append(x0); acc["V1"][1].append(y)
                # V3: moveout as V0, constant loss over the window
                far3 = atten(np.full(N, a0c + axc)) @ conv(moveout(A, tx))
                y = slope(near3, far3, a, b, ta, tb)
                if y is not None:
                    acc["V3"][0].append(np.pi * axc); acc["V3"][1].append(y)
                # V2 / V2n: one event per window
                m1 = moveout(r1, tx); cm1 = conv(m1)
                for v, op, nr in (("V2", atten, near2), ("V2n", atten_amp, near2n)):
                    y = slope(nr, op(np.full(N, a01 + ax1)) @ cm1, a, b, ta, tb)
                    if y is not None:
                        acc[v][0].append(np.pi * ax1); acc[v][1].append(y)
        for v in VARS:
            xs, ys = map(np.array, acc[v]); out[v].append(float(np.log(np.sum(xs ** 2) / np.sum(xs * ys))))
    return fp, out


if __name__ == "__main__":
    t0 = time.monotonic()
    with Pool(len(FPS)) as pool:
        res = dict(pool.map(one, FPS))
    rec = {"fp": {str(int(fp)): {v: dict(per_seed=r[v], bias=float(np.mean(r[v]))) for v in VARS} for fp, r in res.items()},
           "seconds": time.monotonic() - t0}
    p = Path(__file__).resolve().parents[1] / "results" / "diag_resid.json"; p.write_text(json.dumps(rec, indent=1), encoding="utf-8")
    for fp, r in res.items():
        print("fp %2.0f  " % fp + "  ".join("%s %+.3f" % (v, np.mean(r[v])) for v in VARS), flush=True)
