"""Shared code of experiment 10: Q factor from the shallow-to-target spectral decay of the coherent signal.

The experiment-8 generator (semi.py, c395b78) and the experiment-7 library (lib.py, b1d76f8) are imported
unchanged from cache/inputs; semi's paths are redirected to this experiment.

Signal power spectrum of a window: the real part of the adjacent-trace cross-spectrum after removing each
pair's dip lag (C(f), unbiased by noise that is uncorrelated between traces).

Primary estimator SR (`estimate`): K windows of 0.6 s every 0.3 s over 1.2-3.0 s. 0.5 ln C_k(f) is modelled as
a_k + g(f) + 0.5 ln sum_t taper_k(t)^2 exp(-2 pi f t*(t) / c), with t* of the layered model and free window
constants a_k and a free common spectrum g(f) (source); c minimises the rms residual over 10 Hz to the
deepest window's coherent SNR >= 1 limit (cap 80 Hz). Using several windows averages the reflectivity colour
of single windows. Uncertainty: moving-block bootstrap over trace pairs combined with a leave-one-window-out
jackknife (the jackknife carries the window-to-window reflectivity colour); the 90% interval is
c * exp(+-t(0.95, K-1) * hypot(se_boot, se_jack)) on the log scale (t because the jackknife has K-1 degrees
of freedom; chosen on the development seeds, where 1.645 covered 76% and t 85%).

Development variant (`estimate2`, not used for decisions): one shallow (1.2-1.6 s) and one target (2.0-3.0 s)
window.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "cache" / "inputs"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
sys.path.insert(0, str(IN / "novelty_evidence_20261025" / "scripts"))
sys.path.insert(0, str(IN / "field_application_20261028" / "scripts"))
sys.path.insert(0, str(IN / "semisynthetic_volve_20261028" / "scripts"))
import semi  # noqa: E402  (experiment 8, c395b78)
from semi import np, lib, bm, TW, DT  # noqa: E402

semi.ROOT, semi.IN = ROOT, IN            # read_lines reads cache/SNAPSHOT.json; curtain/log files from IN

SHALLOW, TARGET = (1.2, 1.6), (2.0, 3.0)
F_LO, F_CAP = 10.0, 80.0
WL, STEP, T_A, T_B = 0.6, 0.3, 1.2, 3.0
CGRID = np.exp(np.linspace(np.log(0.25), np.log(4.0), 561))     # 0.5% steps
NFFT = 512
T90 = {3: 2.353, 4: 2.132, 5: 2.015, 6: 1.943}                 # t(0.95, df)
MAXLAG = 0.012                                                   # s, dip search per trace pair


def pair_lags(X, maxlag=MAXLAG):
    """Lag (s) of trace j+1 against trace j in a tapered window: cross-correlation peak, parabolic refinement."""
    n = X.shape[0]; F = np.fft.rfft(X, 2 * n, axis=0); cc = np.fft.irfft(F[:, :-1] * np.conj(F[:, 1:]), 2 * n, axis=0)
    m = int(round(maxlag / DT)); cc = np.concatenate([cc[-m - 1:], cc[:m + 2]])    # lags -m-1..m+1
    k = np.argmax(cc[1:-1], axis=0) + 1; j = np.arange(cc.shape[1])
    c0, a, b = cc[k, j], cc[k - 1, j], cc[k + 1, j]; den = a - 2 * c0 + b
    d = np.where(den < 0, 0.5 * (a - b) / np.where(den < 0, den, -1), 0.0)
    # cc[l] = sum x_j(t) x_{j+1}(t - l): peak at l = -(lag of j+1 behind j)
    return -(k - (m + 1) + d) * DT


def spectra(D, t0, t1, align=True, cols=None):
    """f, P (mean power), C (dip-corrected adjacent cross-power), per-pair cross terms (for the bootstrap)."""
    w = semi.win(TW, t0, t1); X = D[w] * lib.taper(w.sum())[:, None]
    F = np.fft.rfft(X, NFFT, axis=0); f = np.fft.rfftfreq(NFFT, DT)
    live = np.abs(D[w]).sum(axis=0) > 0; pair = live[:-1] & live[1:]
    cross = F[:, :-1] * np.conj(F[:, 1:])
    if align:
        cross = cross * np.exp(-2j * np.pi * f[:, None] * pair_lags(X)[None, :])
    return f, np.abs(F) ** 2, np.real(cross), live, pair


def decay_table(cm, Q, f, grid=CGRID, shallow=SHALLOW, target=TARGET):
    """(len(grid), len(f)) predicted 0.5 ln(P_target / P_shallow) for each c."""
    ws, wt = semi.win(TW, *shallow), semi.win(TW, *target)
    ps, pt = lib.taper(ws.sum()) ** 2, lib.taper(wt.sum()) ** 2
    t_s = lib.absorption(dict(cm, Q=np.asarray(Q, float)), TW[ws]); t_t = lib.absorption(dict(cm, Q=np.asarray(Q, float)), TW[wt])
    out = np.empty((len(grid), len(f)))
    for i, c in enumerate(grid):
        a = np.log(ps @ np.exp(-2 * np.pi * np.outer(t_s / c, f))); b = np.log(pt @ np.exp(-2 * np.pi * np.outer(t_t / c, f)))
        out[i] = 0.5 * (b - a)
    return out


def band(f, snr_s, snr_t, fpk):
    """[F_LO, f_hi]: f_hi = lowest of the two windows' contiguous SNR >= 1 runs above the peak, capped."""
    return F_LO, min(lib.snr_band(f, snr_s, fpk), lib.snr_band(f, snr_t, fpk), F_CAP)


def fit(y, table, f, sel):
    """c minimising sum over sel of (y - table_c - const)^2; returns c, rms, index."""
    r = y[None, sel] - table[:, sel]; r = r - r.mean(axis=1, keepdims=True); e = np.sqrt(np.mean(r ** 2, axis=1))
    i = int(np.argmin(e)); return float(CGRID[i]), float(e[i]), i


def estimate2(D, cm, Q, align=True, nboot=0, block=10, seed=0, table=None, shallow=SHALLOW, target=TARGET):
    """SR estimate of the Q factor c with diagnostics; optional moving-block bootstrap over trace pairs."""
    f, Ps, Cs, live_s, pair_s = spectra(D, *shallow, align); _, Pt, Ct, live_t, pair_t = spectra(D, *target, align)
    pair = pair_s & pair_t; live = live_s & live_t

    def curves(ip, il):
        cs = lib.smooth(Cs[:, ip].mean(axis=1), 5); ct = lib.smooth(Ct[:, ip].mean(axis=1), 5)
        ps = lib.smooth(Ps[:, il].mean(axis=1), 5); pt = lib.smooth(Pt[:, il].mean(axis=1), 5)
        snr_s = np.maximum(cs, 0) / np.maximum(ps - np.maximum(cs, 0), 1e-30 * ps.max())
        snr_t = np.maximum(ct, 0) / np.maximum(pt - np.maximum(ct, 0), 1e-30 * pt.max())
        return cs, ct, ps, snr_s, snr_t

    if table is None:
        table = decay_table(cm, Q, f, shallow=shallow, target=target)
    ip = np.where(pair)[0]; il = np.where(live)[0]
    cs, ct, ps, snr_s, snr_t = curves(ip, il)
    fpk = float(f[np.argmax(np.sqrt(ps) * (f >= 5))]); lo, hi = band(f, snr_s, snr_t, fpk)
    sel = (f >= lo) & (f <= hi) & (cs > 0) & (ct > 0)
    y = np.full_like(f, np.nan); y[sel] = 0.5 * np.log(ct[sel] / cs[sel])
    c, rms, _ = fit(y, table, f, sel)
    slope = float(np.polyfit(f[sel], y[sel], 1)[0]) if sel.sum() > 2 else float("nan")
    out = dict(c=c, rms=rms, f_hi=float(hi), n_freq=int(sel.sum()), slope_per_Hz=slope)
    if nboot:
        rng = np.random.default_rng(seed); n = len(ip); cb = []
        for _ in range(nboot):
            st = rng.integers(0, n - block + 1, int(np.ceil(n / block)))
            jj = np.concatenate([ip[s:s + block] for s in st])[:n]
            cs_b = lib.smooth(Cs[:, jj].mean(axis=1), 5); ct_b = lib.smooth(Ct[:, jj].mean(axis=1), 5)
            ok = sel & (cs_b > 0) & (ct_b > 0); yb = np.full_like(f, np.nan); yb[ok] = 0.5 * np.log(ct_b[ok] / cs_b[ok])
            cb.append(fit(yb, table, f, ok)[0])
        cb = np.array(cb); out.update(ci90=[float(np.quantile(cb, 0.05)), float(np.quantile(cb, 0.95))],
                                     boot_sd_log=float(np.std(np.log(cb))))
    out["_curve"] = dict(f=f, y=y, sel=sel)
    return out


# ---------------------------------------------------------------- primary: multi-window

def windows(wl=WL, step=STEP, t0=T_A, t1=T_B):
    return [(float(a), float(a + wl)) for a in np.round(np.arange(t0, t1 - wl + 1e-9, step), 6)]


def mw_table(cm, Q, f, wins, grid=CGRID):
    """(len(grid), K, len(f)) predicted 0.5 ln window power for each c."""
    out = np.empty((len(grid), len(wins), len(f)))
    for k, (a, b) in enumerate(wins):
        w = semi.win(TW, a, b); tp = lib.taper(w.sum()) ** 2; ts = lib.absorption(dict(cm, Q=np.asarray(Q, float)), TW[w])
        for i, c in enumerate(grid):
            out[i, k] = 0.5 * np.log(tp @ np.exp(-2 * np.pi * np.outer(ts / c, f)))
    return out


def mw_fit(L, table, sel, keep=None):
    """c minimising rms of L - table_c after removing window and frequency means (double centring)."""
    keep = np.ones(L.shape[0], bool) if keep is None else keep
    r = L[None][:, keep][:, :, sel] - table[:, keep][:, :, sel]
    r = r - r.mean(axis=1, keepdims=True); r = r - r.mean(axis=2, keepdims=True)
    e = np.sqrt(np.mean(r ** 2, axis=(1, 2))); i = int(np.argmin(e))
    return float(CGRID[i]), float(e[i])


def estimate(D, cm, Q, align=True, nboot=0, block=10, seed=0, table=None, wins=None):
    """SR estimate of the Q factor c (primary); diagnostics and, with nboot, the combined 90% interval."""
    wins = windows() if wins is None else wins
    sp = [spectra(D, a, b, align) for a, b in wins]; f = sp[0][0]
    pair = np.logical_and.reduce([s[4] for s in sp]); live = np.logical_and.reduce([s[3] for s in sp])
    ip, il = np.where(pair)[0], np.where(live)[0]
    if table is None:
        table = mw_table(cm, Q, f, wins)

    def curves(jj):
        return np.array([lib.smooth(s[2][:, jj].mean(axis=1), 5) for s in sp])

    C = curves(ip); P = np.array([lib.smooth(s[1][:, il].mean(axis=1), 5) for s in sp])
    snr = np.maximum(C, 0) / np.maximum(P - np.maximum(C, 0), 1e-30 * P.max(axis=1, keepdims=True))
    fpk = float(f[np.argmax(np.sqrt(np.maximum(C[0], 0)) * (f >= 5))])
    hi = min(lib.snr_band(f, snr[-1], fpk), F_CAP)
    sel = (f >= F_LO) & (f <= hi) & np.all(C > 0, axis=0)
    L = np.zeros_like(C); L[:, sel] = 0.5 * np.log(C[:, sel])
    c, rms = mw_fit(L, table, sel)
    jack = [mw_fit(L, table, sel, np.arange(len(wins)) != k)[0] for k in range(len(wins))]
    lj = np.log(jack); se_jack = float(np.sqrt((len(wins) - 1) / len(wins) * np.sum((lj - lj.mean()) ** 2)))
    out = dict(c=c, rms=rms, f_hi=float(hi), n_freq=int(sel.sum()), n_windows=len(wins), jack=jack,
               se_jack_log=se_jack)
    if nboot:
        rng = np.random.default_rng(seed); n = len(ip); cb = []
        for _ in range(nboot):
            st = rng.integers(0, n - block + 1, int(np.ceil(n / block)))
            Cb = curves(np.concatenate([ip[s:s + block] for s in st])[:n]); ok = sel & np.all(Cb > 0, axis=0)
            Lb = np.zeros_like(Cb); Lb[:, ok] = 0.5 * np.log(Cb[:, ok]); cb.append(mw_fit(Lb, table, ok)[0])
        se_boot = float(np.std(np.log(cb), ddof=1)); se = float(np.hypot(se_boot, se_jack)); k = T90[len(wins) - 1]
        out.update(se_boot_log=se_boot, se_log=se, ci90=[float(c * np.exp(-k * se)), float(c * np.exp(k * se))])
    return out
