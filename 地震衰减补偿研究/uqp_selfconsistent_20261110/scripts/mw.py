"""UQ-MW: multichannel Bayes-Wiener inverse Q with an attenuation-adaptive lateral aperture (exploratory prototype,
NOTE_ZH.md section 3).

Time-frequency cell (t, f), K = 2 m + 1 adjacent traces: d_j = G(u) x_j + e_j, with source-level signal x_j of
power S0 and lateral coherence E[x_i x_j*] = S0 gamma(j - i), noise e_j of power N, uncorrelated between traces.
The linear estimate of the centre trace that minimises the expected (over the attenuation multiplier u) error is

    x_c = E[G*]  sum_j w_j d_j,     w = (E|G|^2 Gamma + (N / S0) I)^-1 gamma_c

  K = 1                    -> E[G*] / (E|G|^2 + N / S0): the single-trace operator of uqp.py;
  gamma = 1 (flat events)  -> the K-fold stack followed by that operator with N / K;
  N / (S0 E|G|^2) -> 0     -> w = unit vector: no lateral mixing where the signal is strong;
  N / (S0 E|G|^2) large    -> w proportional to gamma_c: matched lateral filter.
So the lateral aperture widens with f t*: resolution is traded only where attenuation has pushed the signal under
the noise. The signal power must be predicted from a reference level (uqp.surface_fit): after lateral mixing a
local coherence SNR cannot be measured, because the mixed noise is coherent between neighbours.

Prototype restrictions: one coherence function per (window, frequency) for the whole line (no lateral variation of
dip), real coherence (flat events), geometric decay gamma(d) = q^|d| with q from the fitted lag-2 / lag-1 surfaces.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402
from uqp import np, lib, expectations, UGRID, WIN, HOP, DT, TINY  # noqa: E402


def lag_spectra(D, tw, lags=(1, 2)):
    """Window centres; f; total power P; raw real cross-power at each trace lag (n_lag, n_win, n_f)."""
    cen, P, C = [], [], []
    for a in np.arange(tw[0], tw[-1] - WIN + 1e-9, HOP):
        w = (tw >= a - 1e-9) & (tw <= a + WIN + 1e-9); X = np.fft.rfft(D[w] * lib.taper(w.sum())[:, None], 512, axis=0)
        P.append(lib.smooth(np.mean(np.abs(X) ** 2, axis=1), 5)); cen.append(a + WIN / 2)
        C.append([lib.smooth(np.mean(np.real(X[:, :-d] * np.conj(X[:, d:])), axis=1), 5) for d in lags])
    return np.fft.rfftfreq(512, DT), np.array(cen), np.array(P), np.array(C).transpose(1, 0, 2)


def signal_model(f, cen, P, C, tw, x, dx, u=1.0):
    """Reference-level signal model at attenuation multiplier u: ln S0 (n_win, n_f), lateral coherence q per trace
    lag (geometric), noise power N. Each lag surface is fitted with alpha_w beta(f) <|G_u|^2>_w."""
    i = int(np.argmin(np.abs(UGRID - u))); ug = UGRID[i:i + 1]; tab = uqp.atten_table(tw, x, dx, f, cen, ug); wt = uqp.weights(f, P)
    f1 = uqp.surface_fit(C[0], wt, tab, ugrid=ug); f2 = uqp.surface_fit(C[1], wt, tab, ugrid=ug)
    q = np.clip(np.exp(uqp.log_s0(f2) - uqp.log_s0(f1)), 0.05, 1.0)             # gamma(2) / gamma(1)
    ls0 = uqp.log_s0(f1) - np.log(q)                                            # S0 = S0 gamma(1) / q
    # noise = total minus the locally measured coherent power extrapolated to lag 0 (not minus the model: the local
    # excess of signal over the smooth model must not be counted as noise)
    N = np.clip(P - np.maximum(C[0], 0) / q, 0.02 * P, P)
    return ls0, q, N


def rows(f, cen, arr, tw, f_op, log=True):
    a = np.log(arr) if log else arr; a = np.array([np.interp(f_op, f, r) for r in a])
    a = np.stack([np.interp(tw, cen, a[:, k]) for k in range(len(f_op))], axis=1)
    return np.exp(a) if log else a


def weights_rows(q, rho_c, m):
    """w = (Gamma + rho_c I)^-1 gamma_c for every (t, f): q, rho_c (n, nf) -> (n, nf, 2 m + 1). rho_c = N / (S0 E|G|^2)."""
    k = 2 * m + 1; d = np.abs(np.arange(k)[:, None] - np.arange(k)[None, :])
    G = q[..., None, None] ** d + rho_c[..., None, None] * np.eye(k); g = q[..., None] ** np.abs(np.arange(k) - m)
    return np.linalg.solve(G, g[..., None])[..., 0]


def uqmw(D, tw, x_c, dx_c, sig, ls0, q, N, f, cen, f_op, fref, s, m):
    """x_c = E[G*] sum_j w_j d_j with the lognormal expectation around x_c (relative spread sig on dx_c)."""
    lf = np.log(np.maximum(f_op, np.finfo(float).tiny) / float(fref)); A, Cinv, _ = expectations(x_c, dx_c, sig, f_op, lf, "lognormal")
    rho = uqp.rho_rows(f, cen, N, ls0, tw, f_op); qr = np.clip(rows(f, cen, q, tw, f_op), 0.05, 1.0)
    W = weights_rows(qr, (s + rho) * Cinv, m)                                    # rho_c = (s + N / S0) / E|G|^2
    nt = D.shape[1]; pad = np.pad(D, ((0, 0), (m, m)), mode="reflect") if m else D; out = np.zeros_like(D)
    for j in range(2 * m + 1):
        out += uqp.apply_rows(A * W[..., j], pad[:, j:j + nt])                    # A = E[G*] / E|G|^2
    return out
