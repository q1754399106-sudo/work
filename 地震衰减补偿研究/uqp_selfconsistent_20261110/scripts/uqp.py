"""UQ-P: self-consistent Bayes-Wiener inverse Q (exploratory prototype, NOTE_ZH.md).

One attenuation posterior is used in the three places where attenuation enters a Wiener inverse-Q operator:
  (i)   amplitude and (ii) dispersion of the operator:  E_q[G*] / E_q|G|^2;
  (iii) the signal power that the noise term is measured against: predicted from a reference level by the
        attenuation law, not measured cell by cell.

    H(t, f) = E_q[G*] / (E_q|G|^2 + s + N(t, f) / S0(t, f))

  N(t, f)  incoherent power of the data (P - C of lib.coherence_snr), measurable everywhere;
  S0(t, f) de-attenuated coherent power alpha_w beta(f): a window constant times a common spectrum, from a
           weighted least-squares fit of the raw adjacent-trace cross-power over the whole band;
  q        lognormal posterior of the attenuation multiplier u (true t* = u x assumed t* past the anchor):
           prior ln u ~ N(0, sig0^2) around the supplied Q, likelihood from the same surface fit
           C_w(f) ~ alpha_w beta(f) <|G_u|^2>_w (window jackknife standard error).

Limits: flat likelihood -> prior only (UQ-W's expectation with a reference-level noise term);
        sharp likelihood -> plug-in Wiener inverse Q at the estimated attenuation.
UQ-W (uqw.py) measures the signal power locally instead: rho = s + E|G|^2 / snr(t, f).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "qunknown_comp_20261109" / "scripts"))
from uqw import np, lib, semi, expectations, DT  # noqa: E402

WIN, HOP = 0.4, 0.1
UGRID = np.exp(np.linspace(np.log(0.25), np.log(4.0), 281))      # 1% steps
SIG0 = 0.447                                                     # prior: Q known within a factor of 2
BAND, TINY = (8.0, 80.0), 1e-300


def win_spectra(D, tw):
    """Window centres; f; raw coherent power C (adjacent-trace cross-power, may be negative), incoherent power
    N = P - max(C, 0) and total power P. Same windows, taper and 5-point smoothing as lib.coherence_snr."""
    cen, C, P = [], [], []
    for a in np.arange(tw[0], tw[-1] - WIN + 1e-9, HOP):
        w = (tw >= a - 1e-9) & (tw <= a + WIN + 1e-9); X = np.fft.rfft(D[w] * lib.taper(w.sum())[:, None], 512, axis=0)
        live = np.abs(D[w]).sum(axis=0) > 0; pair = live[:-1] & live[1:]
        P.append(lib.smooth(np.mean(np.abs(X[:, live]) ** 2, axis=1), 5))
        C.append(lib.smooth(np.mean(np.real(X[:, :-1] * np.conj(X[:, 1:]))[:, pair], axis=1), 5)); cen.append(a + WIN / 2)
    C, P = np.array(C), np.array(P)
    return np.fft.rfftfreq(512, DT), np.array(cen), C, P - np.maximum(C, 0), P


def win_atten(tw, x, f, cen):
    """ln of the taper^2-weighted mean of exp(-2 pi f x(t)) over each window: (n_win, n_f)."""
    out = np.empty((len(cen), len(f)))
    for k, c in enumerate(cen):
        w = (tw >= c - WIN / 2 - 1e-9) & (tw <= c + WIN / 2 + 1e-9); tp = lib.taper(w.sum()) ** 2
        l = -2 * np.pi * np.outer(x[w], f); mx = l.max(axis=0)
        out[k] = mx + np.log((tp[:, None] * np.exp(l - mx)).sum(axis=0) / tp.sum())
    return out


def atten_table(tw, x_abs, dx, f, cen, ugrid=UGRID):
    """win_atten for t* = x_abs + dx (u - 1), every u of the grid: (n_u, n_win, n_f)."""
    return np.array([win_atten(tw, x_abs + dx * (u - 1), f, cen) for u in ugrid])


def weights(f, P):
    """Inverse-variance weights of the raw cross-power (variance proportional to P^2) inside the band."""
    return np.where((f >= BAND[0]) & (f <= BAND[1]), 1.0, 0.0)[None, :] / np.maximum(P, TINY) ** 2


def _pos(v):
    return np.maximum(v, 1e-12 * np.abs(v).max() + TINY)


def surface_fit(C, wt, table, ugrid=UGRID, g=None, keep=None, n_iter=40):
    """Weighted least squares of the raw coherent power, C_w(f) ~ alpha_w beta(f) exp(table_u), over u (grid,
    parabolic refinement in ln u). Linear domain and every cell of the band: noise-dominated cells enter with
    their variance and without the bias of taking logs of selected positive estimates.
    g: known ln beta ((n_f,) or (n_win, n_f)) or None (free common spectrum, alternating scale fits).
    keep: window mask (jackknife). Returns ln S0 = ln alpha_w + ln beta(f) at the best u."""
    wt = wt if keep is None else wt * keep[:, None]
    nu, nw, nf = table.shape; cost = np.empty(nu); par = []
    for i in range(nu):
        M = np.exp(table[i]) if g is None else np.exp(table[i] + (g if np.ndim(g) == 2 else g[None, :]))
        al = np.ones(nw)
        be = _pos((wt * C * M).sum(axis=0) / np.maximum((wt * M * M).sum(axis=0), TINY)) if g is None else np.ones(nf)
        for _ in range(n_iter if g is None else 1):
            B = M * be[None, :]; al = _pos((wt * C * B).sum(axis=1) / np.maximum((wt * B * B).sum(axis=1), TINY))
            if g is None:
                A = M * al[:, None]; be = _pos((wt * C * A).sum(axis=0) / np.maximum((wt * A * A).sum(axis=0), TINY))
        cost[i] = (wt * (C - al[:, None] * be[None, :] * M) ** 2).sum(); par.append((al, be))
    i = int(np.argmin(cost)); lu = np.log(ugrid); u_hat = float(ugrid[i])
    if 0 < i < nu - 1:
        den = cost[i - 1] - 2 * cost[i] + cost[i + 1]
        if den > 0:
            u_hat = float(np.exp(lu[i] + 0.5 * (cost[i - 1] - cost[i + 1]) / den * (lu[1] - lu[0])))
    al, be = par[i]; lb = np.log(be)[None, :] if g is None else (g if np.ndim(g) == 2 else g[None, :])
    S = al[:, None] * np.exp(lb + table[i])                                    # predicted coherent power
    rel = (wt * S ** 2).sum(axis=1)                                            # window reliability: sum of (S / P)^2
    return dict(u=u_hat, ls0=np.log(al)[:, None] + lb, rel=rel, cost=float(cost[i] / max((wt * C ** 2).sum(), TINY)),
                edge=bool(nu > 1 and i in (0, nu - 1)))


def posterior(fit, C, wt, table, g=None, sig0=SIG0):
    """Lognormal posterior of u: prior N(0, sig0^2) on ln u, likelihood N(ln u_fit, se^2), se from the
    leave-one-window-out jackknife of the surface fit. Returns (u_post, sig_post, se); the prior if the fit
    sits on the grid edge."""
    if fit["edge"]:
        return 1.0, sig0, float("inf")
    lj = []
    for k in range(len(C)):
        keep = np.ones(len(C), bool); keep[k] = False                       # 4% grid: enough for a standard error
        lj.append(np.log(surface_fit(C, wt, table[::4], ugrid=UGRID[::4], g=g, keep=keep, n_iter=20)["u"]))
    n = len(lj); lj = np.array(lj); se = max(float(np.sqrt((n - 1) / n * np.sum((lj - lj.mean()) ** 2))), 1e-3)
    prec = 1 / sig0 ** 2 + 1 / se ** 2
    return float(np.exp((np.log(fit["u"]) / se ** 2) / prec)), float(prec ** -0.5), se


def log_s0(fit, rel_min=3.0):
    """ln S0 (n_win, n_f); windows whose fitted signal is not reliable take the level of the nearest reliable one."""
    ls0 = fit["ls0"].copy(); idx = np.where(fit["rel"] >= rel_min)[0]
    if len(idx) == 0:
        return ls0
    lvl = ls0.mean(axis=1)
    for k in range(len(ls0)):
        if fit["rel"][k] < rel_min:
            ls0[k] += lvl[idx[np.argmin(np.abs(idx - k))]] - lvl[k]
    return ls0


def rho_rows(f, cen, N, ls0, tw, f_op, n_floor=1e-6):
    """N(t, f) / S0(t, f) at every output row and operator frequency (log-linear interpolation)."""
    lr = np.log(np.maximum(N, n_floor * N.max())) - ls0
    lr = np.array([np.interp(f_op, f, r) for r in lr])
    return np.exp(np.stack([np.interp(tw, cen, lr[:, k]) for k in range(len(f_op))], axis=1))


def apply_rows(H, D):
    n = H.shape[0]; nfft = 2 * (H.shape[1] - 1); lag = np.arange(n); M = np.empty((n, n))
    for i in range(n):
        h = H[i].copy(); h[0] = h[0].real; h[-1] = h[-1].real
        M[i] = np.fft.irfft(h, nfft)[(i - lag) % nfft]
    return M @ D


def operator(x_c, dx_c, sig, rho, f_op, fref, s):
    """H = E[G*] / (E|G|^2 + s + rho): lognormal expectation around x_c with relative spread sig on dx_c."""
    lf = np.log(np.maximum(f_op, np.finfo(float).tiny) / float(fref))
    A, Cinv, _ = expectations(x_c, dx_c, sig, f_op, lf, "lognormal")
    return A / (1 + (s + rho) * Cinv)
