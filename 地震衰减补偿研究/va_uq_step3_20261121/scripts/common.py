"""Shared helpers for va_uq_step3_20261121 (PLAN_ZH.md). NumPy only.

Operator tables are ports of uq_rtm_20261117/scripts/uq_tab.m, gain_s.m and
va_uq_step1_20261121/scripts/va_tab.m (same quadrature, same formulas)."""
import numpy as np

SIG = 0.447


def gh41():
    """probabilists' Gauss-Hermite nodes/weights, 41 points (= numpy hermegauss(41), uq_tab.m gh41)."""
    x, w = np.polynomial.hermite_e.hermegauss(41)
    return x, w


def gain_s(g_db):
    target = 10 ** (g_db / 20)
    lo, hi = -40.0, 0.0
    for _ in range(200):
        s = np.exp((lo + hi) / 2)
        b = -s + np.sqrt(s * s + s)
        g = (b + s) / (b * b + s)
        if g > target:
            lo = (lo + hi) / 2
        else:
            hi = (lo + hi) / 2
    return np.exp((lo + hi) / 2)


def tab(f, x, sig, fv):
    """A = E[G*]/E|G|^2, C = 1/E|G|^2, P = E[exp(-ikx)] at scalar f, array x (two-way t*, s).
    fv scalar: phase reference known (uq_tab / va_tab with fixed f_v); fv=(f1,f2): ln f_v uniform (va_tab)."""
    x = np.atleast_1d(np.asarray(x, float))
    if np.isscalar(fv):
        fb, lam = float(fv), 0.0
    else:
        fb, lam = np.sqrt(fv[0] * fv[1]), np.log(fv[1] / fv[0])
    k = 2 * f * np.log(f / fb)
    if sig == 0:
        e, lw = np.zeros(1), np.zeros(1)
    else:
        e, w = gh41()
        lw = np.log(w / w.sum())
    xs = x[None, :] * np.exp(sig * e)[:, None]
    y = f * xs * lam
    sc = np.ones_like(y)
    nz = y != 0
    sc[nz] = np.sin(y[nz]) / y[nz]
    l2 = -2 * np.pi * f * xs + lw[:, None]
    mx = l2.max(0)
    L2 = mx + np.log(np.exp(l2 - mx).sum(0))
    A = (np.exp((-np.pi * f - 1j * k) * xs + lw[:, None] - L2) * sc).sum(0)
    C = np.exp(-L2)
    P = (np.exp(-1j * k * xs + lw[:, None]) * sc).sum(0)
    return A, C, P


S70 = gain_s(70)


def H(method, f, x, inv_snr=0.0):
    """Operator value of the step-1 methods (va_mig.m va_tables, UQW1 form gamma 1, s = gain_s(70)).
    x is the t* the method acts on (absolute t*/m, or anchored Delta t*/m for VA, VAfix)."""
    spec = {'PH100': (0, 100.0, True), 'UQW1_100': (SIG, 100.0, False), 'UQW1_25': (SIG, 25.0, False),
            'VA': (SIG, (15.0, 35.0), False), 'PH25': (0, 25.0, True), 'VAfix': (SIG, 25.0, False)}[method]
    A, C, P = tab(f, x, spec[0], spec[1])
    if spec[2]:
        return P
    return (A + S70 * P * C + P * inv_snr) / (1 + S70 * C + inv_snr)


def wls_slope(f, y, w=None, quad=False):
    """Weighted LS of y = c0 + c1 f (+ c2 f^2); returns c1 and its formal SE (if w are 1/var)."""
    f = np.asarray(f, float)
    X = np.c_[np.ones_like(f), f] if not quad else np.c_[np.ones_like(f), f, f * f]
    W = np.ones_like(f) if w is None else np.asarray(w, float)
    XtW = X.T * W
    cov = np.linalg.inv(XtW @ X)
    beta = cov @ (XtW @ y)
    return beta[1], np.sqrt(cov[1, 1]), beta


def sine_tapers(n, k):
    """Riedel-Sidorenko sine tapers (unit energy), shape (k, n)."""
    t = np.arange(1, n + 1)
    return np.array([np.sqrt(2 / (n + 1)) * np.sin(np.pi * (j + 1) * t / (n + 1)) for j in range(k)])
