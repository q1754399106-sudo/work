"""UQ-MW with dip steering (exploratory development of mw.py for dipping events, NOTE_ZH.md section 9).

Neighbours are aligned to the centre trace before the multichannel weights are applied:
  1. local dip p(t, x): lag of trace x + 1 against trace x in 0.4 s windows every 0.1 s (cross-correlation peak,
     sq.pair_lags), median over DIP_MED traces and 3 windows, linear in time between window centres;
  2. neighbour x + j is read at t + (sum of p over the traces between x and x + j);
  3. lag spectra, signal model and weights as mw.py, with the aligned neighbours.
Coherence is still one function per (window, frequency) for the whole line.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp, mw  # noqa: E402
from uqp import np, lib, UGRID, WIN, HOP, DT  # noqa: E402

MAXLAG, DIP_MED = 0.012, 9


def pair_lags(X, maxlag=MAXLAG):
    """sq.pair_lags (experiment 10), copied: lag (s) of trace j + 1 against trace j in a tapered window."""
    n = X.shape[0]; F = np.fft.rfft(X, 2 * n, axis=0); cc = np.fft.irfft(F[:, :-1] * np.conj(F[:, 1:]), 2 * n, axis=0)
    m = int(round(maxlag / DT)); cc = np.concatenate([cc[-m - 1:], cc[:m + 2]])
    k = np.argmax(cc[1:-1], axis=0) + 1; j = np.arange(cc.shape[1])
    c0, a, b = cc[k, j], cc[k - 1, j], cc[k + 1, j]; den = a - 2 * c0 + b
    d = np.where(den < 0, 0.5 * (a - b) / np.where(den < 0, den, -1), 0.0)
    return -(k - (m + 1) + d) * DT


def dip_field(D, tw):
    """p (n_t, n_x - 1): lag of trace x + 1 against trace x at every sample."""
    cen, P = [], []
    for a in np.arange(tw[0], tw[-1] - WIN + 1e-9, HOP):
        w = (tw >= a - 1e-9) & (tw <= a + WIN + 1e-9); P.append(pair_lags(D[w] * lib.taper(w.sum())[:, None])); cen.append(a + WIN / 2)
    P = np.array(P); h = DIP_MED // 2; pad = np.pad(P, ((1, 1), (h, h)), mode="edge")
    S = np.median(np.stack([pad[i:i + P.shape[0], j:j + P.shape[1]] for i in range(3) for j in range(DIP_MED)]), axis=0)
    return np.stack([np.interp(tw, cen, S[:, k]) for k in range(S.shape[1])], axis=1)


def aligned(D, tw, p, j):
    """Trace x + j read at t + cumulative dip lag from x to x + j (edge traces: mirrored index, zero extra lag)."""
    if j == 0:
        return D
    nt, nx = D.shape; cp = np.concatenate([np.zeros((nt, 1)), np.cumsum(p, axis=1)], axis=1); out = np.empty_like(D)
    for x in range(nx):
        y = x + j; y = -y if y < 0 else (2 * (nx - 1) - y if y > nx - 1 else y)
        out[:, x] = np.interp(tw + (cp[:, y] - cp[:, x]), tw, D[:, y])
    return out


def lag_spectra(D, A1, A2, tw):
    """As mw.lag_spectra with the aligned lag-1 and lag-2 neighbours A1, A2."""
    cen, P, C = [], [], []
    for a in np.arange(tw[0], tw[-1] - WIN + 1e-9, HOP):
        w = (tw >= a - 1e-9) & (tw <= a + WIN + 1e-9); tp = lib.taper(w.sum())[:, None]
        X, X1, X2 = (np.fft.rfft(Z[w] * tp, 512, axis=0) for Z in (D, A1, A2))
        P.append(lib.smooth(np.mean(np.abs(X) ** 2, axis=1), 5)); cen.append(a + WIN / 2)
        C.append([lib.smooth(np.mean(np.real(X * np.conj(Y)), axis=1), 5) for Y in (X1, X2)])
    return np.fft.rfftfreq(512, DT), np.array(cen), np.array(P), np.array(C).transpose(1, 0, 2)


class Gather:
    """Dip field, aligned neighbours and lag spectra of one gather (independent of the assumed Q)."""

    def __init__(self, D, tw, m_max, steer=True):
        self.D, self.tw = D, tw; p = dip_field(D, tw) if steer else np.zeros((D.shape[0], D.shape[1] - 1)); self.p = p
        self.nb = {j: aligned(D, tw, p, j) for j in range(-m_max, m_max + 1)}
        self.f, self.cen, self.P, self.C = lag_spectra(D, self.nb[1], self.nb[2], tw)

    def model(self, x, dx, u=1.0):
        return mw.signal_model(self.f, self.cen, self.P, self.C, self.tw, x, dx, u)

    def apply(self, x_c, dx_c, sig, model, f_op, fref, s, m):
        ls0, q, N = model; lf = np.log(np.maximum(f_op, np.finfo(float).tiny) / float(fref))
        A, Cinv, _ = uqp.expectations(x_c, dx_c, sig, f_op, lf, "lognormal")
        rho = uqp.rho_rows(self.f, self.cen, N, ls0, self.tw, f_op); qr = np.clip(mw.rows(self.f, self.cen, q, self.tw, f_op), 0.05, 1.0)
        W = mw.weights_rows(qr, (s + rho) * Cinv, m); out = np.zeros_like(self.D)
        for j in range(-m, m + 1):
            out += uqp.apply_rows(A * W[..., j + m], self.nb[j])
        return out
