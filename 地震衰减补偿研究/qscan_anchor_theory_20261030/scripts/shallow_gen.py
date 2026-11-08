"""Semi-synthetic generator and field reader with an explicit observation start t0 (PLAN_WSHALLOW_ZH.md).

The experiment-8/10 recipe (semi.py c395b78, gen.py of experiment 10) uses the module-level axes TW (1.2-3.0 s)
and TF (1.0-3.2 s) in read_lines, reflectivity, tstar_true, model and band_power. Those functions are copied
here line for line with the axes as arguments; nothing else changes. With t0 = 1.2 they must reproduce
gen.prepare / gen.signal exactly (checked by --check). The modelling grid always starts 0.2 s above t0.
"""
import sys, json
from env15 import np, lib, semi, gen

DT, FINE = semi.DT, semi.FINE


def axes(t0, t1=3.0):
    tw = np.round(np.arange(int(round(t0 / DT)), int(round(t1 / DT)) + 1) * DT, 6)
    tf = np.round(np.arange(int(round((t0 - 0.2) / FINE)), 3201) * FINE, 6)
    return tw, tf


def read_lines(lines, t0, t1=3.0):
    """semi.read_lines with the window t0-t1 (same replacement rule for incomplete or dead lines)."""
    snap = json.loads((semi.ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))
    d = np.memmap(snap["cube"], dtype=np.uint8, mode="r"); ns = 1126; tb = 240 + 4 * ns; n = (len(d) - 3600) // tb
    hdr = np.asarray(d[3600:3600 + n * tb]).reshape(n, tb)
    il = hdr[:, 188:192].copy().view(">i4").ravel(); xl = hdr[:, 192:196].copy().view(">i4").ravel()
    w = semi.win(np.arange(ns) * DT, t0, t1); want = np.arange(semi.XL[0], semi.XL[1] + 1, 2)
    out, used = {}, set(lines)
    for L in lines:
        while True:
            sel = np.where((il == L) & (xl >= semi.XL[0]) & (xl <= semi.XL[1]))[0]; sel = sel[np.argsort(xl[sel])]
            if np.array_equal(xl[sel], want):
                X = np.stack([semi.ibm(hdr[j, 240:].copy().view(">u4")) for j in sel], axis=1)[w]
                # dead-trace test on the standard 1.2-3.0 s part, so that the line choice is the same as before
                if np.all(np.abs(X[np.arange(len(X)) >= int(round((1.2 - t0) / DT))]).sum(axis=0) > 0):
                    break
            L2 = L + 2
            while L2 in used:
                L2 += 2
            print(f"inline {L} replaced by {L2}", flush=True); used.add(L2); L = L2
        out[L] = X
    return out


def read_inline(L, t0, t1=3.0):
    """lib13.read_inline with the window t0-t1."""
    from field_kappa2 import L13
    h, x, y, il, xl = L13.cube(); w = semi.win(np.arange(1126) * DT, t0, t1)
    D = np.zeros((int(w.sum()), (L13.XL1 - L13.XL0) // L13.XSTEP + 1))
    for j in np.where((il == L) & (xl >= L13.XL0) & (xl <= L13.XL1) & ((xl - L13.XL0) % L13.XSTEP == 0))[0]:
        D[:, (xl[j] - L13.XL0) // L13.XSTEP] = semi.ibm(h[j, 240:].copy().view(">u4"))[w]
    return D


def reflectivity(tau, rng, tw, tf):
    t_log, rc = semi.log_reflectivity()

    def boot(n):
        s = rng.integers(0, len(rc), int(np.ceil(n / semi.BLOCK)))
        return np.concatenate([np.take(rc, np.arange(k, k + semi.BLOCK), mode="wrap") for k in s])[:n]

    nf, nx = len(tf), tau.shape[1]
    r0 = boot(nf); i0 = int(round((t_log[0] - tf[0]) / FINE)); r0[i0:i0 + len(rc)] = rc
    ind = np.stack([boot(nf) for _ in range(nx)], axis=1)
    k = np.arange(-15, 16); g = np.exp(-0.5 * (k / 5.0) ** 2); g /= g.sum()
    ind = np.stack([np.convolve(np.pad(row, 15, mode="reflect"), g, mode="valid") for row in ind])
    ind *= np.sqrt(np.mean(r0 ** 2) / np.mean(ind ** 2))
    base = np.sqrt(1 - semi.ETA ** 2) * r0[:, None] + semi.ETA * ind
    R = np.zeros((nf, nx))
    for j in range(nx):
        pos = (tf + np.interp(tf, tw, tau[:, j]) - tf[0]) / FINE
        i = np.floor(pos).astype(int); fr = pos - i
        for ii, ww in ((i, 1 - fr), (i + 1, fr)):
            ok = (ii >= 0) & (ii < nf); np.add.at(R[:, j], ii[ok], base[ok, j] * ww[ok])
    return R, r0


def tstar_true(cm, c, tf):
    tt = np.arange(0, 3401) * FINE; a = lib.absorption(cm, tt)
    a = np.convolve(np.pad(a, 100, mode="edge"), np.ones(201) / 201, mode="valid")
    return np.interp(tf, tt, a) / c


def source_spectrum(w_in, f, cm, c, tf):
    """semi.source_spectrum (phase 0) with the true t* at 1.4 s taken on tf (identical values)."""
    x = np.zeros(512); w_in = np.asarray(w_in); x[:101] = w_in[100:]; x[-100:] = w_in[:100]
    f4 = np.fft.rfftfreq(512, DT); A = np.abs(np.fft.rfft(x)) * (DT / FINE)
    S = np.interp(f, f4, A, right=0.0)
    dt_star = np.interp(1.4, tf, tstar_true(cm, c, tf)) - lib.absorption(cm, [1.4])[0]
    return S * np.exp(np.pi * f * dt_star) * np.exp(0j * f)


def model(R, tstar, S, f, fref, tw, tf, N=8192):
    aa = np.clip((120 - f) / 20, 0, 1); aa = 0.5 - 0.5 * np.cos(np.pi * aa)
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / fref); b[0] = 0
    off = np.round((tw[:, None] - tf[None, :]) / FINE).astype(int) % N
    M = np.empty((len(tw), len(tf)))
    for k in range(len(tf)):
        H = S * aa * np.exp(b * tstar[k]); H[-1] = H[-1].real
        M[:, k] = np.fft.irfft(H, N)[off[:, k]]
    return M @ R


def band_power(X, tw, a, b, nfft=512):
    w = semi.win(tw, a, b)
    return lib.smooth(np.mean(np.abs(np.fft.rfft(X[w] * lib.taper(w.sum())[:, None], nfft, axis=0)) ** 2, axis=1), 5)


def prepare(inline, noise_lines, t0):
    """gen.prepare on the window t0-3.0 s. Calibration windows (1.2-1.6 s source power, 2.0-3.0 s SNR) unchanged."""
    tw, tf = axes(t0); cm = semi.curtain_cm(); lines = read_lines([inline] + list(noise_lines), t0); used = list(lines)
    D = lines[used[0]]; lags = semi.local_lags(D); tau = semi.structure(lags); wv = lib.wavelets(D, tw, cm)
    noise = [semi.residual_noise(lines[L], semi.local_lags(lines[L])) for L in used[1:]]
    f, snr, P = lib.coherence_snr(D, tw, 1.2, 1.6); Creal = P * snr / (1 + snr)
    fq, snr_real, _ = lib.coherence_snr(D, tw, 2.0, 3.0)
    return dict(cm=cm, D=D, tau=tau, w_in=wv["w_in"], noise=noise, lines=[int(x) for x in used], Creal=Creal,
                band=(f >= 5) & (f <= 80), snr_real=snr_real, bq=(fq >= 10) & (fq <= 80), tw=tw, tf=tf, t0=t0)


def signal(ctx, c, seed):
    """gen.signal (block-bootstrap reflectivity) on the context's axes."""
    cm, tw, tf = ctx["cm"], ctx["tw"], ctx["tf"]; f = np.fft.rfftfreq(8192, FINE); rng = np.random.default_rng(seed)
    R = reflectivity(ctx["tau"], rng, tw, tf)[0]
    S = source_spectrum(ctx["w_in"], f, cm, c, tf)
    x = model(R, tstar_true(cm, c, tf), S, f, float(cm["fref"]), tw, tf)
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / band_power(x, tw, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + L * N0, tw, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in gen.LEVEL_GRID]
    return dict(x=x, level=float(gen.LEVEL_GRID[int(np.argmin(err))]), R=R, scale=k, S=S, f=f)


if __name__ == "__main__" and "--check" in sys.argv:
    # 1) t0 = 1.2 reproduces gen.prepare / gen.signal; 2) field read on 0.8 s, sliced, equals lib13.read_inline
    import time
    t = time.monotonic(); a = gen.prepare(10227, [10149, 10273, 10355]); b = prepare(10227, [10149, 10273, 10355], 1.2)
    print("lines", a["lines"], b["lines"], "tau", np.abs(a["tau"] - b["tau"]).max(), "noise", np.abs(a["noise"][0] - b["noise"][0]).max())
    ga = gen.signal(a, 1.0, 10001); gb = signal(b, 1.0, 10001)
    print("x max|d| %.3e rel %.3e  level %g %g  scale %g %g" % (np.abs(ga["x"] - gb["x"]).max(), np.abs(ga["x"] - gb["x"]).max() / np.abs(ga["x"]).max(),
          ga["level"], gb["level"], ga["scale"], gb["scale"]))
    from field_kappa2 import L13
    d0 = L13.read_inline(10141); d1 = read_inline(10141, 0.8)
    print("field slice max|d|", np.abs(d1[100:] - d0).max(), d1.shape, f"{time.monotonic() - t:.0f}s")
    e = prepare(10227, [10149, 10273, 10355], 0.8)
    print("t0 0.8: lines", e["lines"], "D", e["D"].shape, "tau", e["tau"].shape, "tau 1.2-3.0 vs standard max|d|",
          np.abs(e["tau"][100:] - a["tau"]).max())
    g = signal(e, 1.0, 10025); print("ext signal", g["x"].shape, "level", g["level"], "scale %.4g" % g["scale"], f"{time.monotonic() - t:.0f}s")
