"""Q-compensated Kirchhoff prestack depth migration with interchangeable gain laws (exploratory, synthetic).

Each (shot, receiver, image point) contribution is taken from the trace filtered with H(f, x), where
x = [t*_s + t*_r] / m is the path t* of the migration Q model (true smooth Q times m). Tables: first-arrival
traveltime T (fast marching, smooth c0) and t* = int ds / (c0 Q) (symmetric traveltime perturbation, Fermat).
Forward model of one path: G = exp[(-pi f - i kf) x u], kf = 2 f ln(f / F0), u = true / nominal attenuation.

Gain laws H(f, x):
  NONE  1                                         PH    exp(i kf x)                (phase only)
  LP    exp[(pi f + i kf) x] * cosine low-pass    STAB  exp(i kf x) (b + s) / (b^2 + s), b = exp(-pi f x)
  WIEN  G* / (|G|^2 + rho)                        UQ    E[G*] / (E|G|^2 + rho), ln u ~ N(0, sig^2)
rho = S70 + E|G|^2 / [fold * snr(t, f)]. snr: recorded traces in 0.3 s windows (window power against an
independent noise record of the same level). fold: the gain in snr from stacking, measured as the snr of the
uncompensated image (image of the data against the image of the noise record) over the snr of the traces.
"""
import json, time
import numpy as np
import numba as nb
import skfmm
import scipy.fft as sfft
from common import *

NFINE = 512; DTF = TREC / NFINE; KBIN = 48
S70 = np.exp(-(0.23 * 70 + 1.63)); SIG = 0.447
GH = np.polynomial.hermite_e.hermegauss(41)
KF = 2 * FREQ * np.log(FREQ / F0)
ZONES = dict(all=(slice(19, 132), slice(25, 170)), target=(slice(53, 98), slice(75, 163)))   # image-grid indices


def tables():
    p = CACHE / "tables.npz"
    if p.exists():
        d = np.load(p); return d["T"], d["A"]
    M = model(); x, z = M["x"], M["z"]; X, Z = np.meshgrid(x, z, indexing="ij"); r0 = 2 * H; eps = 0.5
    T = np.empty((len(IXR), len(IXR), NZ), np.float32); A = np.empty_like(T)
    for i, ix in enumerate(IXR):
        phi = np.hypot(X - x[ix], Z - z[ZSR]) - r0; c = M["cs"]; q = M["qis"]
        tt = lambda sp: np.asarray(skfmm.travel_time(phi, sp, dx=H, order=2)) + r0 / sp[ix, ZSR]
        t0 = tt(c); tp = tt(c / (1 + eps * q)); tm = tt(c / (1 - eps * q))
        T[i] = t0[::2]; A[i] = ((tp - tm) / (2 * eps))[::2]
    np.savez(p, T=T, A=A); return T, A


@nb.njit(parallel=True, cache=True)
def mig_shot(bank, Ts, As, T, A, invm, dxk, tmax, xs, xr, xi, zi, img):
    nr, K, N = bank.shape; nxi, nzi = img.shape
    for ix in nb.prange(nxi):
        for iz in range(3, nzi):
            z = zi[iz] - zi[2]; cs = z / np.sqrt(z * z + (xi[ix] - xs) ** 2); acc = 0.0
            for r in range(nr):
                t = Ts[ix, iz] + T[r, ix, iz]
                if t >= tmax:
                    continue
                a = (As[ix, iz] + A[r, ix, iz]) * invm / dxk; k = int(a)
                if k >= K - 1:
                    k = K - 2; a = K - 1.0
                fa = a - k; u = t / DTF; n = int(u); fu = u - n
                v0 = bank[r, k, n] * (1 - fu) + bank[r, k, n + 1] * fu
                v1 = bank[r, k + 1, n] * (1 - fu) + bank[r, k + 1, n + 1] * fu
                acc += (v0 * (1 - fa) + v1 * fa) * cs * z / np.sqrt(z * z + (xi[ix] - xr[r]) ** 2)
            img[ix, iz] += acc


def expect(x, sig):
    """E[G*] / E|G|^2 and E|G|^2 for x (K,), lognormal u; (K, nf). sig = 0: plain G."""
    x = np.asarray(x, float)[:, None]; b = -np.pi * FREQ + 1j * KF
    if sig == 0:
        return np.exp((np.pi * FREQ + 1j * KF) * x), np.exp(-2 * np.pi * FREQ * x)
    e, w = GH; lw = np.log(w / w.sum()); xs = x[:, :, None] * np.exp(sig * e)[None, None, :]
    l2 = -2 * np.pi * FREQ[None, :, None] * xs + lw; mx = l2.max(-1, keepdims=True)
    L2 = mx[..., 0] + np.log(np.exp(l2 - mx).sum(-1))
    return np.exp(b[None, :, None] * xs + lw - L2[..., None]).sum(-1), np.exp(L2)


def law(name, x, par=None, inv_snr=None):
    x = np.asarray(x, float); xc = x[:, None]
    if name == "NONE":
        return np.ones((len(x), len(FREQ)), complex)
    if name == "PH":
        return np.exp(1j * KF * xc)
    if name == "LP":
        tap = 0.5 * (1 + np.cos(np.pi * np.clip((FREQ - par) / 8.0, 0, 1)))
        return np.exp((np.pi * FREQ + 1j * KF) * xc) * tap
    if name == "STAB":
        s = np.exp(-(0.23 * par + 1.63)); b = np.exp(-np.pi * FREQ * xc)
        return np.exp(1j * KF * xc) * (b + s) / (b * b + s)
    Ac, G2 = expect(x, SIG if par is None else par) if name == "UQ" else expect(x, 0.0)
    if inv_snr is None:
        return Ac / (1 + S70 / G2)
    return Ac[None] / (1 + S70 / G2[None] + inv_snr[:, None, :])                 # (nwin, K, nf)


def to_time(Y):
    """(..., nf) spectra at JF -> real traces (..., NFINE + 1)."""
    buf = np.zeros(Y.shape[:-1] + (NFINE,), complex); buf[..., JF] = Y
    y = 2 * sfft.fft(buf, axis=-1, workers=4).real
    return np.concatenate([y, y[..., :1]], -1)


WIN = int(round(0.3 / DTF)); HOP = int(round(0.25 / DTF)); STARTS = np.arange(0, NFINE - WIN + 1, HOP)
CEN = (STARTS + WIN / 2) * DTF                                                  # window centres (s)


def inv_snr_windows(D, Nz):
    """1/snr of the recorded traces, (nwin, nf). D: noisy data (nf, ns, nr); Nz: an independent noise record of
    the same level (stands for a noise estimate from before the first breaks). 0.3 s Hann windows, all traces."""
    hw = np.hanning(WIN); fw = np.fft.rfftfreq(WIN, DTF)
    d = to_time(np.moveaxis(D, 0, -1))[..., :NFINE]; n = to_time(np.moveaxis(Nz, 0, -1))[..., :NFINE]; inv = []
    for a in STARTS:
        P = (abs(np.fft.rfft(d[..., a:a + WIN] * hw, axis=-1)) ** 2).mean((0, 1))
        Pn = (abs(np.fft.rfft(n[..., a:a + WIN] * hw, axis=-1)) ** 2).mean((0, 1))
        inv.append(np.interp(FREQ, fw, Pn / np.maximum(P - Pn, 1e-3 * Pn)))
    return np.array(inv)


def migrate(D, H, m, xk, T, A, phase=None):
    """D (nf, ns, nr); H (K, nf), or (nwin, K, nf) for filters that vary with record time (blended linearly
    between the window centres) -> image (len(IXR), NZ)."""
    M = model(); xi = M["x"][::2].copy(); zi = M["z"].copy(); img = np.zeros((len(IXR), NZ))
    F = 2 * np.pi * FREQ * (PHASE if phase is None else phase); t = np.arange(NFINE + 1) * DTF
    if H.ndim == 3:
        wt = np.array([np.interp(t, CEN, np.eye(len(CEN))[w]) for w in range(len(CEN))])
    for s, ixs in enumerate(IXS):
        d = D[:, s, :].T[:, None, :]
        if H.ndim == 2:
            bank = to_time(d * (H * F)[None])                                    # (nr, K, NFINE + 1)
        else:
            bank = sum(to_time(d * (H[w] * F)[None]) * wt[w] for w in range(len(CEN)))
        mig_shot(bank, T[ixs // 2], A[ixs // 2], T, A, 1.0 / m, xk[1] - xk[0], TREC - 0.06,
                 xi[ixs // 2], xi, xi, zi, img)
    return img


PHASE = np.exp(1j * float(os.environ.get("QMIG_PHASE", "0")))                 # imaging-filter phase, set by check_phase


def score(I, R):
    out = {}
    for k, (sx, sz) in ZONES.items():
        a = I[sx, sz].ravel(); r = R[sx, sz].ravel()
        out[k] = dict(corr=float(a @ r / np.sqrt((a @ a) * (r @ r) + 1e-300)), E=float(1 - ((a - r) ** 2).sum() / (r @ r)))
    return out
