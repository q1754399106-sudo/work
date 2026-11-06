"""Python port of the operator of qcdn_core.m / time_absorp_wavelet_matrix.m (Wang et al. 2021), so the repo lasso solver can run on it.

A(fw, Q) = W(fw) G(fw, Q), W the Ricker matrix, G the absorption matrix. The port is checked against MATLAB by dev_lasso.py.
"""
import numpy as np


def ricker_matrix(n, dt, fw):
    t = (np.arange(n) - (n - 1) / 2) * dt
    w = (1 - 2 * (np.pi * fw * t) ** 2) * np.exp(-(np.pi * fw * t) ** 2); w = w / np.abs(w).max()
    nn = (n + 1) // 2; M = np.zeros((n, n))
    for i in range(1, nn + 1):
        M[: nn + i - 1, i - 1] = w[nn - i:]
    for i in range(1, nn):
        M[i:, i + nn - 1] = w[: n - i]
    return M


def absorp_matrix(n, dt, fw, Q):
    """n odd, as in the experiment (1001)."""
    df = 1 / (n * dt); nn = (n - 1) // 2 + 1; f_ref = 10 * fw
    f = df * np.arange(1, nn)[:, None]; tk = (np.arange(n) * dt)[None, :]
    H = np.zeros((n, n), complex)
    H[1:nn] = np.exp(-np.pi * f * tk / Q) * np.exp(1j * 2 * np.log(f / f_ref) * f * tk / Q)
    H[nn:] = np.conj(H[1:nn])[::-1]
    T = np.real(np.fft.ifft(H, axis=0)); ss = np.vstack([T[nn:], T[:nn]])
    TW = np.zeros((n, n))
    for i in range(1, nn + 1):
        TW[: nn + i - 1, i - 1] = ss[nn - i:, i - 1]
    for i in range(1, nn):
        TW[i:, i + nn - 1] = ss[: n - i, i + nn - 1]
    return TW


def operator(n, dt, fw, Q):
    W = ricker_matrix(n, dt, fw)
    return W @ absorp_matrix(n, dt, fw, Q), W
