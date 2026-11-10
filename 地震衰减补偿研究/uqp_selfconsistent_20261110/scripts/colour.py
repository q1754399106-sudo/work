"""Depth-stationary colour model for attenuation estimation from stacked data (exploratory, NOTE_ZH.md).

Spectral-ratio estimators compare time windows and assume that the reflectivity spectrum (colour) is the same in
every window. Layering is laid down in depth: if its statistics are stationary in depth and velocity grows with
depth, the reflectivity series in two-way time is compressed with depth and its spectrum is a scaled copy,
R_w(f) = Phi(f / v_w). The drift mimics (or hides) attenuation. Here the coherent power is modelled as

    C_w(f) ~ alpha_w  beta(f)  Phi(f / v_w)  <|G_u|^2>_w(f)

alpha_w: window constants; beta(f): common spectrum in frequency (source, processing); Phi(k): common spectrum in
the scaled frequency k = f / v_w (reflectivity in depth), a free positive function (hat functions in ln k);
v_w: interval velocity of the window from the public velocity model. Phi = const is the usual assumption.
The attenuation multiplier u is identified by the part of t*(t) that is not an affine function of 1 / v(t).

surface_fit3 mirrors uqp.surface_fit (weighted least squares of the raw cross-power, grid over u).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqp import np, lib, UGRID, BAND, TINY, WIN, _pos  # noqa: E402


def window_velocity(tw, v_of_t, cen):
    """1 / (taper^2-weighted mean of 1 / v) over each window."""
    out = []
    for c in cen:
        w = (tw >= c - WIN / 2 - 1e-9) & (tw <= c + WIN / 2 + 1e-9); tp = lib.taper(w.sum()) ** 2
        out.append(tp.sum() / (tp / v_of_t[w]).sum())
    return np.array(out)


def layered_velocity(cm, tw):
    """Interval velocity at two-way time tw for the public layered model (cm['vp'], cm['h'])."""
    v = np.asarray(cm["vp"], float); h = np.asarray(cm["h"], float); tb = 2 * np.cumsum(h / v[:len(h)])
    return v[np.searchsorted(tb, tw, side="right").clip(0, len(v) - 1)]


def scale_basis(f, v_w, n_knots):
    """Hat-function index and weight of ln(f / v_w) on a uniform grid over the band: (n_win, n_f) each."""
    fb = (f >= BAND[0]) & (f <= BAND[1]); y = np.log(np.maximum(f, 1e-6))[None, :] - np.log(v_w)[:, None]
    lo, hi = y[:, fb].min(), y[:, fb].max(); pos = np.clip((y - lo) / (hi - lo) * (n_knots - 1), 0, n_knots - 1 - 1e-9)
    i0 = np.floor(pos).astype(int); return i0, pos - i0


def surface_fit3(C, wt, table, f, v_w, ugrid=UGRID, keep=None, n_knots=12, n_iter=60, smooth=1e-3):
    """As uqp.surface_fit with the extra common factor Phi(f / v_w). Returns u, ln S0 = ln(alpha beta Phi), Phi knots."""
    wt = wt if keep is None else wt * keep[:, None]
    nu, nw, nf = table.shape; i0, w1 = scale_basis(f, v_w, n_knots); cost = np.empty(nu); par = []
    D2 = np.diff(np.eye(n_knots), 2, axis=0); R = D2.T @ D2
    for i in range(nu):
        M = np.exp(table[i]); al = np.ones(nw); c = np.ones(n_knots)
        be = _pos((wt * C * M).sum(axis=0) / np.maximum((wt * M * M).sum(axis=0), TINY))
        for _ in range(n_iter):
            ph = c[i0] * (1 - w1) + c[i0 + 1] * w1
            B = M * ph * be[None, :]; al = _pos((wt * C * B).sum(axis=1) / np.maximum((wt * B * B).sum(axis=1), TINY))
            A = M * ph * al[:, None]; be = _pos((wt * C * A).sum(axis=0) / np.maximum((wt * A * A).sum(axis=0), TINY))
            B = M * al[:, None] * be[None, :]; h0, h1 = B * (1 - w1), B * w1; N_ = np.zeros((n_knots, n_knots)); r = np.zeros(n_knots)
            np.add.at(N_, (i0, i0), wt * h0 * h0); np.add.at(N_, (i0 + 1, i0 + 1), wt * h1 * h1)
            np.add.at(N_, (i0, i0 + 1), wt * h0 * h1); np.add.at(N_, (i0 + 1, i0), wt * h0 * h1)
            np.add.at(r, i0, wt * C * h0); np.add.at(r, i0 + 1, wt * C * h1)
            c = _pos(np.linalg.solve(N_ + smooth * np.trace(N_) / n_knots * R + 1e-12 * np.trace(N_) * np.eye(n_knots), r))
            c = c / np.exp(np.mean(np.log(c)))                                # scale goes to alpha
        ph = c[i0] * (1 - w1) + c[i0 + 1] * w1
        B = M * ph * be[None, :]; al = _pos((wt * C * B).sum(axis=1) / np.maximum((wt * B * B).sum(axis=1), TINY))
        cost[i] = (wt * (C - al[:, None] * be[None, :] * ph * M) ** 2).sum(); par.append((al, be, c))
    i = int(np.argmin(cost)); lu = np.log(ugrid); u_hat = float(ugrid[i])
    if 0 < i < nu - 1:
        den = cost[i - 1] - 2 * cost[i] + cost[i + 1]
        if den > 0:
            u_hat = float(np.exp(lu[i] + 0.5 * (cost[i - 1] - cost[i + 1]) / den * (lu[1] - lu[0])))
    al, be, c = par[i]; ph = c[i0] * (1 - w1) + c[i0 + 1] * w1; ls0 = np.log(al)[:, None] + np.log(be)[None, :] + np.log(ph)
    S = np.exp(ls0 + table[i])
    return dict(u=u_hat, ls0=ls0, rel=(wt * S ** 2).sum(axis=1), phi=c, cost=float(cost[i] / max((wt * C ** 2).sum(), TINY)),
                curve=(cost / max((wt * C ** 2).sum(), TINY)).tolist(), edge=bool(nu > 1 and i in (0, nu - 1)))
