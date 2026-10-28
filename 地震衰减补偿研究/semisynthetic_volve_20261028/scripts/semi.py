"""Shared code of experiment 8 (PLAN_ZH.md): semi-synthetic Volve data with known reflectivity and Q.

Signal: F-11 T2 log reflectivity (1 ms) laid along the structure of inline 10141, attenuated with a smooth
Q profile and the source wavelet estimated from the real line, modelled on a 1 ms grid and decimated.
Noise: real noise from other inlines (lateral second difference after dip alignment) or Gaussian noise
with the same time-variant spectrum. The experiment-7 code (lib.py) and the solver are imported unchanged.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "cache" / "inputs"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
sys.path.insert(0, str(IN / "novelty_evidence_20261025" / "scripts"))
sys.path.insert(0, str(IN / "field_application_20261028" / "scripts"))
import lib  # noqa: E402  (experiment 7, b1d76f8; imports the solver snapshot above)
from lib import np, bm  # noqa: E402
from scipy.io import loadmat  # noqa: E402

DT, FINE = 0.004, 0.001
T0, T1 = 1.2, 3.0
TW = np.round(np.arange(int(round(T0 / DT)), int(round(T1 / DT)) + 1) * DT, 6)
TF = np.round(np.arange(1000, 3201) * FINE, 6)            # 1 ms modelling grid, 1.0-3.2 s
INLINE, XL = 10141, (2100, 2500)
NOISE_LINES = [10041, 10061, 10081, 10101, 10121, 10161, 10181, 10201, 10221, 10241]
C_TRUE = [0.64, 0.8, 1.0, 1.58]
ETA, BLOCK, SEED = 0.3, 40, 20261028
CENTRE = 100                                             # "well" column (xline 2300)
BANDS = [(10, 25), (25, 40), (40, 55), (55, 70), (70, 85)]
WIN_LOG, WIN_TARGET = (2.30, 2.80), (2.0, 3.0)


def win(t, a, b):
    return (t >= a - 1e-9) & (t <= b + 1e-9)


def curtain_cm():
    c = loadmat(IN / "field_vsp_volve_f11_20261026/inputs/posthoc_curtain.mat", simplify_cells=True)
    cm = dict(c["cm"]); cm["Q"] = bm.vec(cm["Q"]); cm["solver_offset"] = 0; cm["t"] = TW.copy(); cm["tw"] = TW.copy()
    return cm


def log_reflectivity():
    r = loadmat(IN / "field_vsp_volve_f11_20261026/analysis/posthoc_references.mat", simplify_cells=True)
    # rc(2) is NaN in the archived reference (first 1 ms bin before the first log sample); it is set to 0.
    return np.round(bm.vec(r["t1ms"]), 6), np.nan_to_num(bm.vec(r["rc"]))


# ---------------------------------------------------------------- cube

def ibm(u):
    u = u.astype(np.uint64); s = 1 - 2 * ((u >> 31) & 1).astype(float)
    e = ((u >> 24) & 127).astype(float) - 64; f = (u & 0xFFFFFF).astype(float) / 16777216
    return s * f * 16.0 ** e


def read_lines(lines):
    """{inline: (451, 201) window 1.2-3.0 s}; a line that is incomplete or has a dead trace is replaced by
    the next odd inline (+2) not already used (PLAN_ZH.md)."""
    snap = json.loads((ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))
    d = np.memmap(snap["cube"], dtype=np.uint8, mode="r"); ns = 1126; tb = 240 + 4 * ns; n = (len(d) - 3600) // tb
    hdr = np.asarray(d[3600:3600 + n * tb]).reshape(n, tb)
    il = hdr[:, 188:192].copy().view(">i4").ravel(); xl = hdr[:, 192:196].copy().view(">i4").ravel()
    w = win(np.arange(ns) * DT, T0, T1); want = np.arange(XL[0], XL[1] + 1, 2)
    out, used = {}, set(lines)
    for L in lines:
        while True:
            sel = np.where((il == L) & (xl >= XL[0]) & (xl <= XL[1]))[0]; sel = sel[np.argsort(xl[sel])]
            if np.array_equal(xl[sel], want):
                X = np.stack([ibm(hdr[j, 240:].copy().view(">u4")) for j in sel], axis=1)[w]
                if np.all(np.abs(X).sum(axis=0) > 0):
                    break
            L2 = L + 2
            while L2 in used:
                L2 += 2
            print(f"inline {L} replaced by {L2}", flush=True); used.add(L2); L = L2
        out[L] = X
    return out


# ---------------------------------------------------------------- structure and real noise

def local_lags(D, wl=50, step=25, maxlag=3):
    """Lag (samples) of trace j+1 against trace j at every sample: tapered 200 ms windows every 100 ms,
    cross-correlation peak within +-3 samples with parabolic refinement, 5-trace lateral smoothing."""
    n, nx = D.shape; tp = lib.taper(wl); cen, rows = [], []
    for s in range(0, n - wl + 1, step):
        a = D[s:s + wl] * tp[:, None]
        cc = np.stack([np.sum(a[max(0, -k):wl - max(0, k), :-1] * a[max(0, k):wl - max(0, -k), 1:], axis=0)
                       for k in range(-maxlag - 1, maxlag + 2)])          # lags -4..4
        k = np.argmax(cc[1:-1], axis=0) + 1                               # peak within -3..3
        j = np.arange(nx - 1); c0, cm1, cp1 = cc[k, j], cc[k - 1, j], cc[k + 1, j]
        den = cm1 - 2 * c0 + cp1
        d = np.where(den < 0, 0.5 * (cm1 - cp1) / np.where(den < 0, den, -1), 0.0)
        rows.append(np.clip(k - (maxlag + 1) + d, -maxlag, maxlag)); cen.append(s + (wl - 1) / 2)
    L = np.array(rows)
    L = np.stack([lib.smooth(np.pad(r, 2, mode="edge"), 5)[2:-2] for r in L])
    return np.stack([np.interp(np.arange(n), cen, L[:, j]) for j in range(nx - 1)], axis=1)


def sinc_shift(x, pos, half=8):
    """x sampled at integer indices, evaluated at fractional positions (Hann-windowed sinc, 2*half taps)."""
    i0 = np.floor(pos).astype(int); fr = pos - i0; out = np.zeros(len(pos))
    for m in range(-half + 1, half + 1):
        k = np.clip(i0 + m, 0, len(x) - 1); u = m - fr
        out += x[k] * np.sinc(u) * (0.5 + 0.5 * np.cos(np.pi * u / half))
    return out


def residual_noise(D, lags):
    """Real noise: (d_j - (aligned d_{j-1} + aligned d_{j+1}) / 2) * sqrt(2/3); edge traces use the
    one-sided difference / sqrt(2)."""
    n, nx = D.shape; t = np.arange(n, dtype=float); E = np.empty_like(D)
    for j in range(nx):
        prev = sinc_shift(D[:, j - 1], t - lags[:, j - 1]) if j > 0 else None
        nxt = sinc_shift(D[:, j + 1], t + lags[:, j]) if j < nx - 1 else None
        if prev is not None and nxt is not None:
            E[:, j] = (D[:, j] - 0.5 * (prev + nxt)) * np.sqrt(2 / 3)
        else:
            E[:, j] = (D[:, j] - (prev if prev is not None else nxt)) / np.sqrt(2)
    return E


def gaussian_noise(E, rng, edges=(1.2, 1.65, 2.1, 2.55, 3.0), nfft=1024):
    """Spatially white Gaussian noise with the time-variant amplitude spectrum and rms of the real noise E:
    one white field filtered by each window's spectrum, scaled to that window's rms, blended linearly
    between window centres."""
    n, nx = E.shape; Z = rng.standard_normal((n, nx)); Fz = np.fft.rfft(Z, nfft, axis=0)
    parts, cen = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        w = win(TW, a, b); A = np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(E[w] * lib.taper(w.sum())[:, None],
                                                                               nfft, axis=0)) ** 2, axis=1), 5))
        F = np.fft.irfft(Fz * A[:, None], nfft, axis=0)[:n]
        F *= np.sqrt(np.mean(E[w] ** 2) / np.mean(F[w] ** 2)); parts.append(F); cen.append((a + b) / 2)
    out = np.zeros_like(E)
    for i, F in enumerate(parts):
        lam = np.interp(TW, cen, np.eye(len(cen))[i]); out += lam[:, None] * F
    return out


def structure(lags):
    """tau(t, x): time of an event of the centre trace in trace x, relative to the centre (seconds)."""
    n, nx = lags.shape[0], lags.shape[1] + 1; tau = np.zeros((n, nx))
    for j in range(CENTRE + 1, nx):
        tau[:, j] = tau[:, j - 1] + lags[:, j - 1]
    for j in range(CENTRE - 1, -1, -1):
        tau[:, j] = tau[:, j + 1] - lags[:, j]
    return tau * DT


# ---------------------------------------------------------------- reflectivity and modelling

def reflectivity(tau, rng):
    """(len(TF), nx) at 1 ms: sqrt(1-eta^2) * r0 + eta * (laterally smoothed independent series), each trace
    displaced by tau. r0 = F-11 T2 log RC at its own times, circular 40 ms block bootstrap of it elsewhere."""
    t_log, rc = log_reflectivity()

    def boot(n):
        s = rng.integers(0, len(rc), int(np.ceil(n / BLOCK)))
        return np.concatenate([np.take(rc, np.arange(k, k + BLOCK), mode="wrap") for k in s])[:n]

    nf, nx = len(TF), tau.shape[1]
    r0 = boot(nf); i0 = int(round((t_log[0] - TF[0]) / FINE)); r0[i0:i0 + len(rc)] = rc
    ind = np.stack([boot(nf) for _ in range(nx)], axis=1)
    k = np.arange(-15, 16); g = np.exp(-0.5 * (k / 5.0) ** 2); g /= g.sum()
    ind = np.stack([np.convolve(np.pad(row, 15, mode="reflect"), g, mode="valid") for row in ind])
    ind *= np.sqrt(np.mean(r0 ** 2) / np.mean(ind ** 2))
    base = np.sqrt(1 - ETA ** 2) * r0[:, None] + ETA * ind
    R = np.zeros((nf, nx))
    for j in range(nx):
        pos = (TF + np.interp(TF, TW, tau[:, j]) - TF[0]) / FINE        # deposit with linear split
        i = np.floor(pos).astype(int); fr = pos - i
        for ii, ww in ((i, 1 - fr), (i + 1, fr)):
            ok = (ii >= 0) & (ii < nf); np.add.at(R[:, j], ii[ok], base[ok, j] * ww[ok])
    return R, r0


def tstar_true(cm, c):
    """Two-way t* at TF: layered model t* smoothed over 200 ms (continuous Q), divided by c (Q_true = c * Q)."""
    tt = np.arange(0, 3401) * FINE; a = lib.absorption(cm, tt)
    a = np.convolve(np.pad(a, 100, mode="edge"), np.ones(201) / 201, mode="valid")
    return np.interp(TF, tt, a) / c


def source_spectrum(w_in, f, phase_deg=0.0, cm=None, c=1.0):
    """Source of scenario c on f: the estimated input wavelet (4 ms, centred, 201 samples; the shallow
    1.2-1.6 s spectrum divided by the model absorption at 1.4 s) re-corrected to the scenario's true
    absorption at 1.4 s, so that every scenario reproduces the real shallow spectrum. Zero beyond 125 Hz,
    scaled by DT/FINE so that the 1 ms model has the amplitude of the 4 ms wavelet."""
    x = np.zeros(512); w_in = np.asarray(w_in); x[:101] = w_in[100:]; x[-100:] = w_in[:100]
    f4 = np.fft.rfftfreq(512, DT); A = np.abs(np.fft.rfft(x)) * (DT / FINE)
    S = np.interp(f, f4, A, right=0.0)
    if cm is not None:
        dt_star = np.interp(1.4, TF, tstar_true(cm, c)) - lib.absorption(cm, [1.4])[0]
        S = S * np.exp(np.pi * f * dt_star)
    return S * np.exp(1j * np.deg2rad(phase_deg) * (f > 0))


def model(R, tstar, S, f, fref, N=8192):
    """Attenuated traces at TW: reflector k at TF[k] carries irfft(S * AA * exp(b * t*_k)); 1 ms modelling,
    100-120 Hz anti-alias taper, sampled at 4 ms. tstar=None gives the unattenuated response."""
    aa = np.clip((120 - f) / 20, 0, 1); aa = 0.5 - 0.5 * np.cos(np.pi * aa)
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / fref); b[0] = 0
    off = np.round((TW[:, None] - TF[None, :]) / FINE).astype(int) % N
    M = np.empty((len(TW), len(TF)))
    for k in range(len(TF)):
        H = S * aa * (np.exp(b * tstar[k]) if tstar is not None else 1); H[-1] = H[-1].real
        M[:, k] = np.fft.irfft(H, N)[off[:, k]]
    return M @ R


def cut(f, fc):
    """5 Hz cosine high cut ending at fc and nothing below 3 Hz (as lib.wavelets)."""
    g = np.ones_like(f); r = (f > fc - 5) & (f < fc)
    g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fc - 5)) / 5)); g[f >= fc] = 0; g[f < 3] = 0
    return g


def band_power(X, a, b, nfft=512):
    w = win(TW, a, b)
    return lib.smooth(np.mean(np.abs(np.fft.rfft(X[w] * lib.taper(w.sum())[:, None], nfft, axis=0)) ** 2, axis=1), 5)


def true_band(sig, noise, nfft=512):
    """Highest frequency of the contiguous run above the signal peak with true SNR >= 1 in 2.0-3.0 s."""
    f = np.fft.rfftfreq(nfft, DT); ps = band_power(sig, *WIN_TARGET); pn = band_power(noise, *WIN_TARGET)
    fpk = float(f[np.argmax(np.sqrt(ps) * (f >= 5))])
    return lib.snr_band(f, ps / np.maximum(pn, 1e-30 * ps.max()), fpk)


# ---------------------------------------------------------------- conventional inverse Q (baseline)

def gain_sigma2(g_db):
    """sigma^2 of the stabilised inverse-Q amplitude (beta + s) / (beta^2 + s) whose maximum gain is g_db."""
    target = 10 ** (g_db / 20); lo, hi = -40.0, 0.0
    for _ in range(200):
        s = np.exp((lo + hi) / 2); b = -s + np.sqrt(s * s + s); g = (b + s) / (b * b + s)
        lo, hi = ((lo + hi) / 2, hi) if g > target else (lo, (lo + hi) / 2)
    return float(np.exp((lo + hi) / 2))


def inverse_q(D, cm, Q, g_db, nfft=1024):
    """Time-variant stabilised inverse-Q filter (amplitude and dispersion), model t* of the layered Q."""
    c = dict(cm); c["Q"] = np.asarray(Q, float); ts = lib.absorption(c, TW); n = len(TW)
    f = np.fft.rfftfreq(nfft, DT); s = gain_sigma2(g_db); M = np.empty((n, n)); lag = np.arange(n)
    for i in range(n):
        beta = np.exp(-np.pi * f * ts[i]); lam = (beta + s) / (beta ** 2 + s)
        ph = np.exp(-2j * f * ts[i] * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])))
        H = lam * ph; H[0] = H[0].real; H[-1] = H[-1].real
        M[i] = np.fft.irfft(H, nfft)[(i - lag) % nfft]
    return M @ D


def lowpass(X, fc, nfft=1024):
    f = np.fft.rfftfreq(nfft, DT)
    return np.fft.irfft(np.fft.rfft(X, nfft, axis=0) * cut(f, fc)[:, None], nfft, axis=0)[:X.shape[0]]


def bandpass(X, a, b, nfft=1024):
    f = np.fft.rfftfreq(nfft, DT); g = np.zeros_like(f)
    up = (f > a - 5) & (f < a); g[up] = 0.5 - 0.5 * np.cos(np.pi * (f[up] - (a - 5)) / 5)
    g[(f >= a) & (f <= b)] = 1; dn = (f > b) & (f < b + 5); g[dn] = 0.5 + 0.5 * np.cos(np.pi * (f[dn] - b) / 5)
    return np.fft.irfft(np.fft.rfft(X, nfft, axis=0) * g[:, None], nfft, axis=0)[:X.shape[0]]


# ---------------------------------------------------------------- scores

def colcorr(A, B):
    a = A - A.mean(axis=0); b = B - B.mean(axis=0)
    return np.sum(a * b, axis=0) / np.sqrt(np.sum(a * a, axis=0) * np.sum(b * b, axis=0) + 1e-300)


def ricker_ref(r_fine, fp):
    """Log-synthetic analogue at the well column: RC (1 ms) * Ricker(fp), sampled at 4 ms (as score_field.m)."""
    k = np.arange(-100, 101) * FINE; z = np.pi * fp * k; ric = (1 - 2 * z ** 2) * np.exp(-z ** 2)
    return np.interp(TW, TF, np.convolve(r_fine, ric, mode="same"))


def scores(out, ideal, src, tie_ref):
    """M1 (2.30-2.80 s) and M1w (2.0-3.0 s): median trace correlation with the ideal target; M2: median
    trace correlation per band against the unattenuated source response in 2.0-3.0 s; T2-like tie."""
    w1, w2 = win(TW, *WIN_LOG), win(TW, *WIN_TARGET)
    s = dict(M1=float(np.median(colcorr(out[w1], ideal[w1]))), M1w=float(np.median(colcorr(out[w2], ideal[w2]))))
    for a, b in BANDS:
        s[f"M2_{a}_{b}"] = float(np.median(colcorr(bandpass(out, a, b)[w2], bandpass(src, a, b)[w2])))
    x = out[:, CENTRE - 2:CENTRE + 3].mean(axis=1)
    s["T2like"], s["T2like_shift_ms"] = lib.tie(x, tie_ref, TW, w1)
    _, s["band05_Hz"], s["band95_Hz"] = lib.spectral(out, out, TW)
    return s
