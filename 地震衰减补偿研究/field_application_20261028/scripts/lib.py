"""Shared code of experiment 7 (PLAN_ZH.md): real-data adaptation of SC / QScan + SC.

Changes against the experiment-6 pipeline (the solver itself is imported unchanged):
  * input wavelet estimated from the data (zero phase, amplitude spectrum of the shallow window with the
    model absorption removed, inside the shallow SNR band) instead of a Ricker wavelet;
  * output wavelet = input wavelet limited to the band where the lateral-coherence SNR of the target
    window is >= 1, so that compensation cannot rebuild frequencies the data do not carry;
  * QScan objective (lasso, as qscan_select.m) re-implemented on this dictionary, wider grid 0.4-2.5,
    with a boundary flag.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAP = ROOT / "cache" / "inputs" / "novelty_evidence_20261025"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(SNAP / "scripts"))
import benchmark as bm  # noqa: E402  (common.py pins BLAS to one thread before NumPy loads)
import numpy as np  # noqa: E402

SEL_SC = (0.125, 10.0)          # SC a, mu_relative (novelty_evidence selection.json)
SEL_LASSO = 0.5                 # QScan lasso weight, as qscan_select.m
QGRID = np.exp(np.linspace(np.log(0.4), np.log(2.5), 15))
DT = 0.004


def taper(n, r=0.2):
    w = np.ones(n)
    m = int(np.floor(r * (n - 1) / 2))
    if m > 0:
        k = np.arange(m + 1)
        w[:m + 1] = 0.5 * (1 - np.cos(np.pi * k / m))
        w[-(m + 1):] = w[:m + 1][::-1]
    return w


def smooth(x, k):
    return np.convolve(x, np.ones(k) / k, mode="same")


def coherence_snr(D, tw, t0, t1, nfft=512):
    """SNR(f) = C/(P-C): C = real mean cross-power of adjacent traces, P = mean power (5-point smoothing)."""
    w = (tw >= t0 - 1e-9) & (tw <= t1 + 1e-9)
    X = np.fft.rfft(D[w] * taper(w.sum())[:, None], nfft, axis=0)
    live = np.abs(D[w]).sum(axis=0) > 0
    P = smooth(np.mean(np.abs(X[:, live]) ** 2, axis=1), 5)
    pair = live[:-1] & live[1:]
    C = smooth(np.mean(np.real(X[:, :-1] * np.conj(X[:, 1:]))[:, pair], axis=1), 5)
    f = np.fft.rfftfreq(nfft, DT)
    snr = np.maximum(C, 0) / np.maximum(P - np.maximum(C, 0), 1e-30 * P.max())
    return f, snr, P


def snr_band(f, snr, fpk):
    """Highest frequency of the contiguous run above fpk where SNR >= 1."""
    i = int(np.argmin(np.abs(f - fpk)))
    if snr[i] < 1:
        return float(f[i])
    while i + 1 < len(f) and snr[i + 1] >= 1:
        i += 1
    return float(f[i])


def absorption(cm, t):
    """Two-way cumulative t* at times t for the zero-offset layered model (as common.dictionary)."""
    q = bm.vec(cm["Q"]); v = bm.vec(cm["vp"]); h = bm.vec(cm["h"])
    bounds = np.r_[0, np.cumsum(2 * h / v[:3]), np.inf]
    tau = np.maximum(0, np.minimum(np.asarray(t)[:, None], bounds[1:]) - bounds[:-1])
    return np.sum(tau / q, axis=1)


def wavelets(D, tw, cm, nfft=512):
    """Zero-phase input and output wavelets (time-domain, length 2*100+1 at 4 ms) and band diagnostics."""
    f, snr_s, P_s = coherence_snr(D, tw, 1.20, 1.60, nfft)
    _, snr_d, _ = coherence_snr(D, tw, 2.00, 3.00, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * absorption(cm, [1.40])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))])
    fs = snr_band(f, snr_s, fpk)
    fd = snr_band(f, snr_d, fpk)

    def cut(fc):
        g = np.ones_like(f)
        r = (f > fc - 5) & (f < fc)
        g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fc - 5)) / 5))
        g[f >= fc] = 0
        g[f < 3] = 0
        return g

    a_in = amp * cut(fs)
    a_out = a_in * cut(min(fd, fs))
    k = np.arange(-100, 101)

    def zp(a):
        w = np.fft.irfft(a, nfft)
        w = np.roll(w, 100)[:201]
        return w / np.max(np.abs(w))

    w_in, w_out = zp(a_in), zp(a_out)
    w_out *= np.max(np.abs(np.fft.rfft(w_in, nfft))) / np.max(np.abs(np.fft.rfft(w_out, nfft)))
    return dict(w_in=w_in, w_out=w_out, lag=k * DT, f=f, snr_shallow=snr_s, snr_target=snr_d,
                f_peak=fpk, fmax_shallow=fs, fmax_target=fd)


def toeplitz_wavelet(w, n):
    """W[i, k] = w(t_i - t_k) for a centred wavelet of length 201."""
    W = np.zeros((n, n))
    for k in range(n):
        lo, hi = max(0, k - 100), min(n, k + 101)
        W[lo:hi, k] = w[lo - k + 100:hi - k + 100]
    return W


def dictionary(cm, w_in, w_out):
    """common.dictionary with the Ricker matrix replaced by the input wavelet (B) and output wavelet (P)."""
    t = bm.vec(cm["t"]); n = len(t)
    a = absorption(cm, t)
    N = int(cm["nfft"]); f = np.arange(N // 2 + 1) / (N * float(cm["dt"]))
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); b[0] = 0
    G = np.empty((n, n))
    for k in range(n):
        H = np.exp(b * a[k]); H[-1] = H[-1].real
        G[:, k] = np.fft.irfft(H * np.exp(-2j * np.pi * np.arange(N // 2 + 1) * k / N), n=N)[:n]
    Win, Wout = toeplitz_wavelet(w_in, n), toeplitz_wavelet(w_out, n)
    B = G @ Win
    scale = np.linalg.norm(B, axis=0)
    return B / scale, Wout / scale, scale


def solve(cm, Q, wv, D, sigma, K, live, a, mu):
    c = dict(cm); c["Q"] = np.asarray(Q, float)
    B, P, scale = dictionary(c, wv["w_in"], wv["w_out"])
    out, info = bm.solve(B, P, D, sigma, scale, K, a, mu, live)
    return out, info, (B, P, scale)


def qscan(cm, Q, wv, D, sigma, K, live):
    """Lasso objective of qscan_select.m on the new dictionary over QGRID; returns factor and diagnostics."""
    n = D.shape[0]
    J = []
    for c in QGRID:
        cc = dict(cm); cc["Q"] = np.asarray(Q, float) * c
        B, P, scale = dictionary(cc, wv["w_in"], wv["w_out"])
        # mu = 0, so P only maps the result; the identity returns the coefficients X themselves.
        X, info = bm.solve(B, np.eye(n), D, sigma, scale, K, SEL_LASSO, 0.0, live)
        res = D - B @ X
        pen = SEL_LASSO * np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]
        s2 = np.where(live, sigma ** 2, np.inf)
        J.append(float(np.sum((0.5 * np.sum(res ** 2, axis=0) + np.sum(pen * np.abs(X), axis=0))[live] / s2[live])))
    best = int(np.argmin(J))
    return float(QGRID[best]), dict(grid=QGRID.tolist(), objective=J, best=best,
                                    at_boundary=best in (0, len(QGRID) - 1))


def path_trace(A, idx):
    idx = np.asarray(idx, int) - 1
    return np.array([A[n, idx[n] - 2:idx[n] + 3].mean() for n in range(A.shape[0])])


def zpfilter(x, dt, fhi):
    x = np.asarray(x, float).ravel(); n = len(x); N = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.abs(np.fft.fftfreq(N, dt))
    H = 1 / np.sqrt(1 + (6 / np.maximum(f, 1e-300)) ** 6) / np.sqrt(1 + (f / fhi) ** 24); H[0] = 0
    return np.real(np.fft.ifft(np.fft.fft(x, N) * H))[:n]


def tie(x, ref, tw, w):
    """Max Pearson correlation over +-20 ms shifts (1 ms, linear interpolation), as score_field.m."""
    best, lag = -np.inf, np.nan
    tf, rr = tw[w], ref[w]
    for s in np.arange(-20, 21) / 1000:
        y = np.interp(tf + s, tw, x, left=np.nan, right=np.nan)
        ok = np.isfinite(y) & np.isfinite(rr)
        a, b = y[ok] - y[ok].mean(), rr[ok] - rr[ok].mean()
        c = a @ b / np.sqrt((a @ a) * (b @ b))
        if c > best:
            best, lag = c, s
    return float(best), float(lag * 1000)


def tie_boot(x, xin, ref, tw, w, lag, lag_in, nboot=2000, block=25, seed=20261028):
    """Moving-block bootstrap (blocks of 100 ms) of corr(x) - corr(input) at the fixed best lags."""
    tf, rr = tw[w], ref[w]
    y = np.interp(tf + lag / 1000, tw, x); yi = np.interp(tf + lag_in / 1000, tw, xin)
    n = len(rr); rng = np.random.default_rng(seed); nb = int(np.ceil(n / block)); d = []
    for _ in range(nboot):
        st = rng.integers(0, n - block + 1, nb)
        ii = np.concatenate([np.arange(s, s + block) for s in st])[:n]
        c = lambda u: np.corrcoef(u[ii], rr[ii])[0, 1]  # noqa: E731
        d.append(c(y) - c(yi))
    return [float(v) for v in np.percentile(d, [5, 95])]


def spectral(A, I, tw):
    w = (tw >= 2.0 - 1e-9) & (tw <= 3.0 + 1e-9); n = w.sum(); f = np.arange(n) / (n * DT); b = (f >= 60) & (f <= 100)
    FA = np.abs(np.fft.fft(A[w], axis=0)) ** 2; FI = np.abs(np.fft.fft(I[w], axis=0)) ** 2
    r = FA[b].sum() / FI[b].sum()
    h = f <= 125; p = np.cumsum(FA[h].sum(axis=1)); p /= p[-1]
    return float(r), float(f[h][np.argmax(p >= .05)]), float(f[h][np.argmax(p >= .95)])
