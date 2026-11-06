"""Frequency-saturated offset-consistent damping: generalised Tikhonov (G'G + lam L'L) s = G'd, one step.
L is column-varying and zero-phase with symbol exp[-pi min(f, fb(t, x)) a_x(t, x)]. Below fb the damping follows the offset path at every frequency,
so rho(f, x) = rho(f, 0) exactly in the stationary limit; above fb it stops following, so far-offset high-frequency noise is not boosted.
fb(t) = ln((1/rho_c - 1)/lam) / (2 pi a0(t)): the frequency at which the zero-offset transfer rho(f, a0) falls to rho_c. No wavelet spectrum needed."""
import numpy as np
from gate1 import C, n, NX

N = C.nfft; FF = np.arange(N // 2 + 1) / (N * C.dt); E = np.fft.rfft(np.eye(n), N, axis=0)


def fb_of(a0, lam, rho_c):
    return np.minimum(FF[-1], np.maximum(np.log((1 / rho_c - 1) / lam), 0.0) / (2 * np.pi * np.maximum(a0, 1e-6)))


def Lmat(ax, fb):
    """ax, fb: arrays of length n (one offset)."""
    H = np.exp(-np.pi * np.minimum(FF[:, None], fb[None, :]) * ax[None, :]); return np.fft.irfft(H * E, N, axis=0)[:n]


def fb_floor(fb, ax, wn):
    """Never let the damping at (t, x) fall below wn: exp(-2 pi fb a_x) >= wn."""
    return np.minimum(fb, -np.log(np.maximum(wn, 1e-300)) / (2 * np.pi * np.maximum(ax, 1e-9)))
