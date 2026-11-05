"""Along-horizon spectral estimate of the Q factor (HSR prototype, STAGE0_ZH.md), and the colour-trend generator."""
from env0 import np, lib, semi, sq, fq, TW, DT

CGRID = np.exp(np.linspace(np.log(0.25), np.log(4.0), 281))      # 1% steps
MIN_PAIRS = 20
NGROUPS = 3


def structure(D):
    return semi.structure(semi.local_lags(D))


def window_shifts(tau, a, b):
    """Median tau per trace inside window (a, b) of the centre trace, and the mask of traces whose shifted
    window stays inside the section."""
    w = semi.win(TW, a, b); tb = np.median(tau[w], axis=0)
    valid = (a + tb >= TW[0] - 1e-9) & (b + tb <= TW[-1] + 1e-9)
    return tb, valid


def spread(tb, valid):
    """S_k of STAGE0_ZH.md."""
    if valid.mean() < 0.5:
        return 0.0
    q = np.quantile(tb[valid], [0.1, 0.9]); return float(q[1] - q[0])


def flatten(D, tb):
    """Each trace advanced by its window shift: output sample t holds D(t + tb)."""
    n = D.shape[0]; F = np.fft.rfft(D, 2 * n, axis=0); f = np.fft.rfftfreq(2 * n, DT)
    return np.fft.irfft(F * np.exp(2j * np.pi * f[:, None] * tb[None, :]), 2 * n, axis=0)[:n]


def window_groups(D, tau, a, b, Q_cm):
    """Group spectra of one horizon window: list of (tb_mean, smoothed C(f)), and n pairs per group."""
    cm, Q = Q_cm
    tb, valid = window_shifts(tau, a, b)
    if valid.sum() < 2:
        return []
    edges = np.quantile(tb[valid], np.linspace(0, 1, NGROUPS + 1)[1:-1])
    grp = np.searchsorted(edges, tb); grp[~valid] = -1
    X = flatten(D, tb)
    w = semi.win(TW, a, b); W = X[w] * lib.taper(w.sum())[:, None]
    F = np.fft.rfft(W, sq.NFFT, axis=0); f = np.fft.rfftfreq(sq.NFFT, DT)
    cross = F[:, :-1] * np.conj(F[:, 1:]) * np.exp(-2j * np.pi * f[:, None] * sq.pair_lags(W)[None, :])
    out = []
    for g in range(NGROUPS):
        pair = (grp[:-1] == g) & (grp[1:] == g)
        if pair.sum() < MIN_PAIRS:
            continue
        out.append(dict(tb=float(np.mean(tb[:-1][pair])), C=lib.smooth(np.real(cross[:, pair]).mean(axis=1), 5),
                        pairs=int(pair.sum())))
    return out


def model_table(cm, Q, a, b, tb, f):
    """(len(CGRID), len(f)) predicted 0.5 ln window power for the window shifted by tb."""
    w = semi.win(TW, a, b); tp = lib.taper(w.sum()) ** 2
    ts = lib.absorption(dict(cm, Q=np.asarray(Q, float)), TW[w] + tb)
    return np.stack([0.5 * np.log(tp @ np.exp(-2 * np.pi * np.outer(ts / c, f))) for c in CGRID])


def estimate(D, cm, Q, f_hi, tau=None):
    """HSR estimate; f_hi comes from sq.estimate on the same section (same band as SR)."""
    tau = structure(D) if tau is None else tau
    f = np.fft.rfftfreq(sq.NFFT, DT); band = (f >= sq.F_LO) & (f <= f_hi); fb = f[band]
    sse = np.zeros(len(CGRID)); count = 0; used = []
    for a, b in sq.windows():
        grp = window_groups(D, tau, a, b, (cm, Q))
        if len(grp) < 2:
            continue
        C = np.array([g["C"][band] for g in grp]); ok = np.all(C > 0, axis=0)
        if ok.sum() < 5:
            continue
        L = 0.5 * np.log(C[:, ok])
        T = np.stack([model_table(cm, Q, a, b, g["tb"], fb[ok]) for g in grp], axis=1)   # (nc, ng, nf)
        r = L[None] - T; r = r - r.mean(axis=2, keepdims=True); r = r - r.mean(axis=1, keepdims=True)
        sse += np.sum(r ** 2, axis=(1, 2)); count += r.shape[1] * r.shape[2]
        used.append(dict(win=[a, b], groups=len(grp), n_freq=int(ok.sum()), tb=[g["tb"] for g in grp],
                         pairs=[g["pairs"] for g in grp]))
    if not used:
        return dict(c=float("nan"), used=[])
    e = np.sqrt(sse / count); i = int(np.argmin(e))
    # Curvature of the misfit at the minimum (in log c), a rough measure of how well c is constrained
    j0, j1 = max(i - 5, 0), min(i + 5, len(CGRID) - 1); lc = np.log(CGRID)
    curv = float(np.polyfit(lc[j0:j1 + 1], sse[j0:j1 + 1] / count, 2)[0])
    edge = bool(i == 0 or i == len(CGRID) - 1)
    return dict(c=float(CGRID[i]), rms=float(e[i]), curvature=curv, at_grid_edge=edge, used=used)


def tilt(base, s, f0=35.0, tref=2.1, frame=200, hop=100):
    """Colour trend: amplitude spectrum x exp(s (t_g - tref) (f - f0)), framed along geological time on the
    1 ms grid TF, then rescaled to the original variance."""
    if s == 0:
        return base
    TF = fq.TF; nf = base.shape[0]; h = 0.5 - 0.5 * np.cos(2 * np.pi * (np.arange(frame) + 0.5) / frame)
    nfft = 2 * frame; f = np.fft.rfftfreq(nfft, fq.FINE)
    out = np.zeros_like(base); norm = np.zeros(nf)
    for st in range(-hop, nf, hop):
        i0, i1 = max(st, 0), min(st + frame, nf)
        if i1 - i0 < 2:
            continue
        hh = h[i0 - st:i1 - st]; tc = TF[np.clip(st + frame // 2, 0, nf - 1)]
        G = np.exp(s * (tc - tref) * (f - f0))[:, None]
        x = np.fft.irfft(np.fft.rfft(base[i0:i1] * hh[:, None], nfft, axis=0) * G, nfft, axis=0)[:i1 - i0]
        out[i0:i1] += x * hh[:, None]; norm[i0:i1] += hh ** 2
    out /= np.maximum(norm, 1e-6)[:, None]
    return out * np.sqrt(np.mean(base ** 2) / np.mean(out ** 2))
