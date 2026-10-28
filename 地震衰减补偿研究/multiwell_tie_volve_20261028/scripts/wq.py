"""Shared code of experiment 9: Q factor from the well log (WR), with the reflectivity colour measured, not fitted.

The experiment-10 estimator code (sq.py, gen.py at 346b562), the experiment-8 generator (semi.py, c395b78) and the
experiment-7 library (lib.py, b1d76f8) are imported unchanged from cache/inputs; their paths are redirected here.

Well log -> reflectivity in two-way time (`well`):
  * AC (us/ft) with a 5-point running median; gaps up to GAP_MAX are bridged (linear in MD), longer gaps split
    the log and the longest piece is used;
  * MD -> TVD (and E, N) from the pick file, linear beyond the end picks; vertical two-way sonic time, stretched piecewise
    linearly to the picks that carry a TWT inside the log; outside the outermost such picks the sonic slope is kept;
  * impedance from Gardner (rho ~ v^0.25), RC = diff(ln Z)/2 on a 1 ms grid; position (E, N) along the hole
    from the pick coordinates.

WR estimator (`wr`): windows of 0.6 s spread evenly (step <= 0.3 s) over the log's time span inside 1.2-3.0 s.
For window k, C_k(f) is the dip-corrected adjacent-trace cross-power (sq.spectra) of the 2*HALF+1 traces of the
inline nearest to the hole at the window centre, centred on the nearest trace; R_k(f) is the power of the log RC
in the same window with the same taper. Both are smoothed by SM bins and
    0.5 ln C_k - 0.5 ln R_k = a_k + g(f) + 0.5 ln sum_t taper_k(t)^2 exp(-2 pi f t*(t) / c)
is fitted exactly as SR (sq.mw_fit, double centring), over 10 Hz to the lowest window's coherent SNR >= 1 limit
(cap 80 Hz). Uncertainty as SR: moving-block bootstrap over trace pairs (per window) combined with a
leave-one-window-out jackknife, interval c * exp(+-t(0.95, K-1) * INFL * se).
"""
import sys, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "cache" / "inputs"
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts",
          "semisynthetic_volve_20261028/scripts", "specslope_q_volve_20261028/scripts"):
    sys.path.insert(0, str(IN / p))
import sq, gen  # noqa: E402  (experiment 10, 346b562)
from sq import np, lib, semi, bm, TW, DT  # noqa: E402

sq.ROOT, sq.IN = ROOT, IN
semi.ROOT, semi.IN = ROOT, IN            # read_lines reads cache/SNAPSHOT.json; curtain/log files from IN

E9 = IN / "multiwell_tie_volve_20261028"
GAP_MAX = 40.0                           # m
WL, STEP_MAX = 0.6, 0.3
HALF = 100                               # traces each side of the hole (section step 2 xlines; 100 = whole inline)
XL0, XSTEP = semi.XL[0], 2               # read_lines: xlines 2100-2500 every 2
SM = 5                                   # extra smoothing (bins of 0.49 Hz); development choice, see PLAN_ZH.md
INFL = 1.0                               # interval inflation; development choice, see PLAN_ZH.md
NBOOT, BLOCK = 200, 5
FINE = 0.001


# ---------------------------------------------------------------- wells

def picks():
    """{'15/9-F-12': [(surface, MD, TVD, TWT s or nan, E, N), ...]} for every row with MD and TVD."""
    out, spans = {}, None

    def fl(x):
        try:
            return float(x)
        except ValueError:
            return np.nan
    for line in (E9 / "raw" / "picks" / "Well_picks_Volve_v1.dat").read_text(encoding="latin-1").splitlines():
        if line.strip().startswith("-----"):
            spans = [(m.start(), m.end()) for m in re.finditer(r"-+", line)]; continue
        if spans is None or not line.strip() or line.startswith("Well NO") or line.lstrip().startswith("Well name"):
            continue
        cols = [line[a:b + 1].strip() if i < len(spans) - 1 else line[a:].strip() for i, (a, b) in enumerate(spans)]
        if len(cols) < 12:
            continue
        md, tvd, twt, e, n = fl(cols[4]), fl(cols[5]), fl(cols[7]) / 1000, fl(cols[10]), fl(cols[11])
        if np.isfinite(md) and np.isfinite(tvd):
            out.setdefault(cols[0].replace("NO ", ""), []).append((cols[1], md, tvd, twt, e, n))
    return out


def running_median(x, k=5):
    h = k // 2; xp = np.pad(x, h, mode="edge")
    return np.median(np.stack([xp[i:i + len(x)] for i in range(k)]), axis=0)


def interp_x(x, xp, fp):
    """np.interp with linear extrapolation from the end segments."""
    y = np.interp(x, xp, fp)
    lo, hi = x < xp[0], x > xp[-1]
    y[lo] = fp[0] + (x[lo] - xp[0]) * (fp[1] - fp[0]) / (xp[1] - xp[0])
    y[hi] = fp[-1] + (x[hi] - xp[-1]) * (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
    return y


def well(name, P=None):
    """dict(t: 1 ms TWT grid, rc: Gardner RC (len(t) - 1, between samples), md, E, N on t, t0, t1, anchors)."""
    P = picks() if P is None else P
    A = np.loadtxt(E9 / "inventory" / "logs" / f"15_9-{name}.csv", delimiter=",", skiprows=1)
    md, ac = A[:, 0], A[:, 1]; ok = np.isfinite(ac); md, ac = md[ok], ac[ok]
    cut = np.where(np.diff(md) > GAP_MAX)[0] + 1; pieces = np.split(np.arange(len(md)), cut)
    s = max(pieces, key=lambda p: md[p[-1]] - md[p[0]]); md, ac = md[s], running_median(ac[s])
    pk = sorted(P["15/9-" + name.replace("_", " ")], key=lambda r: r[1])
    pm = np.array([r[1] for r in pk]); pv = np.array([r[2] for r in pk]); pt = np.array([r[3] for r in pk])
    pe = np.array([r[4] for r in pk]); pn = np.array([r[5] for r in pk])
    _, u = np.unique(pm, return_index=True); pm, pv, pt, pe, pn = pm[u], pv[u], pt[u], pe[u], pn[u]
    tvd = interp_x(md, pm, pv); v = 0.3048 / (ac * 1e-6)
    ts = np.r_[0, np.cumsum(2 * np.diff(tvd) / (0.5 * (v[1:] + v[:-1])))]
    a = np.isfinite(pt) & (pm >= md[0]) & (pm <= md[-1]); km, kt = pm[a], pt[a]
    ks = np.interp(km, md, ts); t = ts + np.interp(ts, ks, kt - ks)        # stretch to anchors, sonic slope outside
    tt = np.round(np.arange(np.ceil(t[0] / FINE), np.floor(t[-1] / FINE) + 1) * FINE, 6)
    z = 1.25 * np.log(np.interp(tt, t, v)); mdt = np.interp(tt, t, md)
    g = np.isfinite(pe) & np.isfinite(pn)
    return dict(name=name, t=tt, rc=0.5 * np.diff(z), md=mdt, E=interp_x(mdt, pm[g], pe[g]),
                N=interp_x(mdt, pm[g], pn[g]), t0=float(tt[0]), t1=float(tt[-1]),
                anchors=[[float(x), float(y)] for x, y in zip(km, kt)],
                sonic_minus_pick_ms=[float(1000 * (np.interp(x, md, ts) - np.interp(km[0], md, ts) + kt[0] - y))
                                     for x, y in zip(km, kt)])


def well_windows(w, wl=WL, step=STEP_MAX, lo=1.2, hi=3.0):
    """K windows of wl s spread evenly over the log span inside [lo, hi], start step <= step, on the 4 ms grid."""
    a = max(w["t0"], lo); b = min(w["t1"] - FINE, hi)
    a = TW[np.searchsorted(TW, a - 1e-9)]; b = TW[np.searchsorted(TW, b + 1e-9) - 1]
    if b - a < wl - 1e-9:
        return []
    K = int(np.ceil(round((b - a - wl) / step, 9))) + 1
    st = np.round(np.round(np.linspace(a, b - wl, K) / DT) * DT, 6)
    return [(float(x), float(round(x + wl, 6))) for x in st]


def rc_power(w, a, b, nfft=2048):
    """Power of the log RC in [a, b] (Hann-type taper as the seismic window) on the 512-point 4 ms frequency grid."""
    tm = 0.5 * (w["t"][1:] + w["t"][:-1]); m = (tm >= a - 1e-9) & (tm <= b + 1e-9); r = w["rc"][m]
    X = np.fft.rfft(r * lib.taper(len(r)), nfft); return np.abs(X[: sq.NFFT // 2 + 1]) ** 2


def rc_energy(w, wins, smooth_s=0.04):
    """Local RC energy (running mean of RC^2 over smooth_s) at the TW samples of each window."""
    tm = 0.5 * (w["t"][1:] + w["t"][:-1]); n = int(round(smooth_s / FINE)) | 1
    e = np.convolve(w["rc"] ** 2, np.ones(n) / n, mode="same")
    return [np.interp(TW[semi.win(TW, a, b)], tm, e) for a, b in wins]


def table(cm, Q, f, wins, energy=None, grid=sq.CGRID):
    """sq.mw_table, optionally with the window weights taper^2 multiplied by the local log RC energy."""
    if energy is None:
        return sq.mw_table(cm, Q, f, wins, grid)
    out = np.empty((len(grid), len(wins), len(f)))
    for k, (a, b) in enumerate(wins):
        w = semi.win(TW, a, b); tp = lib.taper(w.sum()) ** 2 * energy[k]
        ts = lib.absorption(dict(cm, Q=np.asarray(Q, float)), TW[w])
        for i, c in enumerate(grid):
            out[i, k] = 0.5 * np.log(tp @ np.exp(-2 * np.pi * np.outer(ts / c, f)))
    return out


def position(w, t):
    return float(np.interp(t, w["t"], w["E"])), float(np.interp(t, w["t"], w["N"]))


# ---------------------------------------------------------------- estimator

def window_spectra(S, a, b):
    """f, mean power (live traces), per-pair cross-power, pair index of one window's section S."""
    f, P, C, live, pair = sq.spectra(S, a, b, True)
    return f, P[:, live].mean(axis=1), C[:, pair]


def wr(secs, R, wins, cm, Q, sm=SM, nboot=0, block=BLOCK, seed=0, table=None, infl=INFL):
    """WR estimate of c from window sections secs[k] (time x traces, TW grid) and log RC powers R[k]."""
    sp = [window_spectra(S, a, b) for S, (a, b) in zip(secs, wins)]; f = sp[0][0]; K = len(wins)
    if table is None:
        table = sq.mw_table(cm, Q, f, wins)
    Rs = np.array([lib.smooth(r, sm) for r in R])

    def curves(idx):
        return np.array([lib.smooth(lib.smooth(s[2][:, j].mean(axis=1), 5), sm) for s, j in zip(sp, idx)])

    full = [np.arange(s[2].shape[1]) for s in sp]
    C5 = np.array([lib.smooth(s[2].mean(axis=1), 5) for s in sp]); P5 = np.array([lib.smooth(s[1], 5) for s in sp])
    snr = np.maximum(C5, 0) / np.maximum(P5 - np.maximum(C5, 0), 1e-30 * P5.max(axis=1, keepdims=True))
    fpk = float(f[np.argmax(np.sqrt(np.maximum(C5[0], 0)) * (f >= 5))])
    hi = min(min(lib.snr_band(f, s, fpk) for s in snr), sq.F_CAP)
    C = curves(full); sel = (f >= sq.F_LO) & (f <= hi) & np.all(C > 0, axis=0) & np.all(Rs > 0, axis=0)

    def logs(Cx, ok):
        L = np.zeros_like(Cx); L[:, ok] = 0.5 * np.log(Cx[:, ok]) - 0.5 * np.log(Rs[:, ok]); return L

    L = logs(C, sel); c, rms = sq.mw_fit(L, table, sel)
    jack = [sq.mw_fit(L, table, sel, np.arange(K) != k)[0] for k in range(K)]
    lj = np.log(jack); se_jack = float(np.sqrt((K - 1) / K * np.sum((lj - lj.mean()) ** 2)))
    out = dict(c=c, rms=rms, f_hi=float(hi), n_freq=int(sel.sum()), n_windows=K, windows=wins, jack=jack,
               se_jack_log=se_jack)
    if nboot:
        rng = np.random.default_rng(seed); cb = []
        for _ in range(nboot):
            idx = []
            for s in sp:
                n = s[2].shape[1]; st = rng.integers(0, n - block + 1, int(np.ceil(n / block)))
                idx.append(np.concatenate([np.arange(x, x + block) for x in st])[:n])
            Cb = curves(idx); ok = sel & np.all(Cb > 0, axis=0); cb.append(sq.mw_fit(logs(Cb, ok), table, ok)[0])
        se_boot = float(np.std(np.log(cb), ddof=1)); se = float(np.hypot(se_boot, se_jack))
        k = sq.T90[K - 1] if K - 1 in sq.T90 else 1.645
        out.update(se_boot_log=se_boot, se_log=se, t90=k,
                   ci90=[float(c * np.exp(-k * infl * se)), float(c * np.exp(k * infl * se))])
    out["_L"], out["_sel"], out["_f"] = L, sel, f
    return out


def pooled(rows, infl=INFL):
    """Inverse-variance pooled log c of several wells; t multiplier with sum(K - 1) degrees of freedom."""
    lc = np.array([np.log(r["c"]) for r in rows]); se = np.array([r["se_log"] for r in rows]); w = 1 / se ** 2
    m = float(np.sum(w * lc) / w.sum()); s = float(1 / np.sqrt(w.sum())); df = int(sum(r["n_windows"] - 1 for r in rows))
    k = {3: 2.353, 4: 2.132, 5: 2.015, 6: 1.943, 7: 1.895, 8: 1.860, 9: 1.833, 10: 1.812}.get(df, 1.645)
    return dict(c=float(np.exp(m)), se_log=s, df=df, t90=k, ci90=[float(np.exp(m - k * infl * s)), float(np.exp(m + k * infl * s))])
