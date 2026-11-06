"""Noise-consistent lower bound for the offset-consistent damping, and lambda selection.
w(t, x) = max( exp(-2 pi fc*(t) a_x), w_n ), where w_n is the smallest w whose closed-form noise gain does not exceed gamma^2 times
the zero-offset noise gain at the same time: N^2(a, w) = sum beta^2 / (beta^2 + lam w)^2, beta = exp(-pi f a)."""
import numpy as np
from gate1 import C, FS, n, NX, DT

F = FS[::8]; LOGW = np.linspace(0, np.log(1e-6), 120)          # 1 down to 1e-6
TSUB = np.arange(0, n, 5)


def noise2(a, lam, w):
    b2 = np.exp(-2 * np.pi * F[None, :] * a[:, None]); return np.sum(b2 / (b2 + lam * w[:, None]) ** 2, axis=1)


def bound(a0, ax, lam, gamma=1.0):
    """w_n on (n, NX): the smallest allowed w. Coarse time grid, linear interpolation in log w."""
    out = np.ones((n, NX)); a0s = a0[TSUB]; N0 = noise2(a0s, lam, np.ones(len(TSUB)))
    for i in range(NX):
        a = a0s + ax[TSUB, i]; lw = np.zeros(len(TSUB))
        for k, t in enumerate(TSUB):
            vals = noise2(np.full(len(LOGW), a[k]), lam, np.exp(LOGW)); ok = vals <= gamma ** 2 * N0[k]
            lw[k] = LOGW[np.where(ok)[0][-1]] if ok.any() else 0.0
        out[:, i] = np.exp(np.interp(np.arange(n), TSUB, lw))
    return out
