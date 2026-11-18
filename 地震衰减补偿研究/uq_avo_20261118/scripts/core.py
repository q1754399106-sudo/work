"""NumPy port of the original split method's model (original_fixed_split_20261022/attenuation_rebuild):
config_model, ray_absorption, attenuation_matrix, make_reference. Same formulas, same default parameters."""
import numpy as np


class C:
    dt = 0.002; nt = 601; t = np.arange(601) * 0.002
    offset = np.arange(0, 1501, 150.0)
    vp = np.array([2100.0, 2250, 2450, 2800]); vs = np.array([550.0, 700, 850, 1000]); rho = np.array([1.8, 2.1, 2.3, 2.4])
    Q = np.array([40.0, 50, 70, 100]); event_t = np.array([0.38, 0.65, 0.90])
    h = np.diff(np.r_[0, event_t]) * vp[:3] / 2
    fpeak = 60.0; fref = 100.0; nfft = 4096


def ray(t0, x, Q=C.Q):
    """a0 (vertical two-way t*), ax (extra t* of the offset ray), p, for times t0 and offset x."""
    t0 = np.atleast_1d(t0); a0 = np.zeros(len(t0)); ax = np.zeros(len(t0)); p = np.zeros(len(t0))
    bound = np.r_[0, np.cumsum(2 * C.h / C.vp[:3]), np.inf]
    for k, T in enumerate(t0):
        tau = np.maximum(0, np.minimum(T, bound[1:]) - bound[:-1]); on = tau > 0
        if not on.any():
            ax[k] = x / (C.vp[0] * Q[0]); continue
        v, q, tt = C.vp[on], Q[on], tau[on]; a0[k] = np.sum(tt / q)
        if x == 0:
            continue
        lo, hi = 0.0, (1 - 1e-12) / v.max()
        for _ in range(70):
            pm = 0.5 * (lo + hi)
            if np.sum(tt * v ** 2 * pm / np.sqrt(1 - (pm * v) ** 2)) < x:
                lo = pm
            else:
                hi = pm
        p[k] = 0.5 * (lo + hi); ax[k] = np.sum(tt / q / np.sqrt(1 - (p[k] * v) ** 2)) - a0[k]
    return a0, ax, p


def atten(a):
    """Column j = response of an impulse at t_j with path loss a[j]."""
    N = C.nfft; f = np.arange(N // 2 + 1) / (N * C.dt); b = -np.pi * f + 2j * f * np.log(np.maximum(f, 1e-300) / C.fref); b[0] = 0
    H = np.exp(np.outer(b, a)); H[-1] = H[-1].real
    E = np.fft.rfft(np.eye(C.nt), N, axis=0)
    return np.fft.irfft(H * E, N, axis=0)[:C.nt]


def ricker(tt, fp=C.fpeak):
    z = np.pi * fp * tt; return (1 - 2 * z ** 2) * np.exp(-z ** 2)


def reference(x):
    """Ideal NMO-aligned trace with Shuey two-term AVO, and the reflection coefficients and (p v)^2 per event."""
    _, _, p = ray(C.event_t, x); s = np.zeros(C.nt); co = np.zeros(3); sin2 = np.zeros(3)
    for j in range(3):
        vp = C.vp[j:j + 2].mean(); vs = C.vs[j:j + 2].mean(); rho = C.rho[j:j + 2].mean()
        dvp = np.diff(C.vp[j:j + 2])[0] / vp; dvs = np.diff(C.vs[j:j + 2])[0] / vs; drho = np.diff(C.rho[j:j + 2])[0] / rho
        A = 0.5 * (dvp + drho); B = 0.5 * dvp - 2 * (vs / vp) ** 2 * (drho + 2 * dvs)
        sin2[j] = (p[j] * C.vp[j]) ** 2; co[j] = A + B * sin2[j]; s += co[j] * ricker(C.t - C.event_t[j])
    return s, co, sin2
