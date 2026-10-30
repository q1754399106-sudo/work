"""New semi-synthetic Volve data with normal faults and a complex, laterally varying Q (DATA_ZH.md).

Same ingredients as experiments 8-14 (F-11 T2 log reflectivity, real structure, real shallow wavelet, real
noise, SNR calibrated to the real line), but on inlines never used before, and with:
  - three normal faults per section (dipping planes, throw growing with depth, missing section, damage zone);
  - Q attached to the faulted geology: the model-Q step moved, random layer perturbations of log Q, a Q
    multiplier per fault block, a lateral trend, and a low-Q gas lens, so true t* is NOT c x model t*;
  - modelling with a t* profile per trace (kernels on a t* grid, linear interpolation between them).
Variant Z is the complex world; its twin Y has the same faults, reflectivity, noise and scale, but the old
recipe's laterally constant Q = c_eff x model Q (no gas), so that later tests can isolate the Q complexity.
"""
import json
from envf import np, lib, semi, gen

TW, TF, FINE, DT = semi.TW, semi.TF, semi.FINE, semi.DT
NX, CENTRE = 201, semi.CENTRE
STRUCT = {"A": 10243, "B": 10291, "C": 10345}
NOISE = {"A": (9987, 10311), "B": (9993, 10327), "C": (10339, 10361)}
GEO_SEED = {"A": 15001, "B": 15002, "C": 15003}
RC_SEED = {(k, r): 15100 + 10 * i + r for i, k in enumerate(STRUCT) for r in (0, 1)}
C0 = (0.75, 1.05, 1.4)
Q_GAS, DTS, HALF = 18.0, 1e-4, 600
FMAX_FIT = (10.0, 60.0)


# ---------------------------------------------------------------- geology

def draw_geology(key, rng):
    """Faults, Q layering, block multipliers, lateral trend and gas lens of one structure (seeded)."""
    while True:
        x0 = np.sort(rng.uniform(25, 175, 3))
        g = rng.uniform(20, 40, 3) * rng.choice([-1, 1], 3)            # traces per second of two-way time
        T = rng.uniform(0.020, 0.060, 3)                                # throw (s) at 3.2 s
        tt = np.linspace(2.2, 2.9, 15)
        well_ok = all(np.min(np.abs(x0[f] + g[f] * (tt - 2.1) - CENTRE)) > 8 for f in range(3))
        if np.min(np.diff(x0)) >= 35 and well_ok:
            break
    faults = [dict(x0=float(x0[f]), g=float(g[f]), throw=float(T[f])) for f in range(3)]
    edges = [1.0]                                                       # random layering of log Q
    while edges[-1] < 3.2:
        edges.append(edges[-1] + float(rng.uniform(0.06, 0.25)))
    vals = np.clip(rng.normal(0, 0.25, len(edges) - 1), -0.5, 0.5)
    lay = dict(edges=[float(e) for e in edges], logq=[float(v) for v in vals], step_shift=float(rng.uniform(-0.15, 0.15)))
    blocks = {code: float(np.clip(np.exp(rng.normal(0, 0.2)), 0.7, 1.4)) for code in range(8)}
    trend = float(rng.uniform(-0.15, 0.15))
    while True:                                                         # gas lens away from the well column
        xg, tgc = rng.uniform(30, 170), rng.uniform(2.05, 2.35)
        W, H = rng.uniform(15, 25), rng.uniform(0.035, 0.06)
        dist = min(abs(f["x0"] + f["g"] * (tgc - 2.1) - xg) for f in faults)
        if abs(xg - CENTRE) >= 30 and dist >= 0.6 * W:
            break
    gas = dict(x=float(xg), t=float(tgc), half_width=float(W), half_thick=float(H), Q=Q_GAS)
    return dict(key=key, faults=faults, layers=lay, blocks=blocks, trend=trend, gas=gas)


def fault_x(f, t):
    return f["x0"] + f["g"] * (t - 2.1)


def throw(f, t):
    return f["throw"] * (0.4 + 0.6 * (t - 1.0) / 2.2)


def hanging(f, t, j):
    """True on the hanging-wall (downthrown) side; the plane dips towards sign(g)."""
    x = fault_x(f, t)
    return (j > x) if f["g"] > 0 else (j < x)


def mapping(geo, tau):
    """Per trace on the 1 ms output grid TF: geological time tg(t, j) (centre-trace frame), block code,
    distance (traces) to the nearest fault plane. Normal faults give upward jumps of tg (missing section)."""
    nf = len(TF); tg = np.empty((nf, NX)); code = np.zeros((nf, NX), int); dist = np.full((nf, NX), np.inf)
    for j in range(NX):
        tj = TF - np.interp(TF, TW, tau[:, j])
        for b, f in enumerate(geo["faults"]):
            h = hanging(f, TF, j); tj = tj - throw(f, TF) * h; code[:, j] += h * (1 << b)
            dist[:, j] = np.minimum(dist[:, j], np.abs(j - fault_x(f, TF)))
        tg[:, j] = np.maximum.accumulate(tj)
    return tg, code, dist


def base_series(rng):
    """sqrt(1-eta^2) r0 + eta (laterally smoothed independent series), on the centre-trace time grid TF
    (semi.reflectivity before its structural deposit)."""
    t_log, rc = semi.log_reflectivity()

    def boot(n):
        s = rng.integers(0, len(rc), int(np.ceil(n / semi.BLOCK)))
        return np.concatenate([np.take(rc, np.arange(k, k + semi.BLOCK), mode="wrap") for k in s])[:n]

    nf = len(TF); r0 = boot(nf); i0 = int(round((t_log[0] - TF[0]) / FINE)); r0[i0:i0 + len(rc)] = rc
    ind = np.stack([boot(nf) for _ in range(NX)], axis=1)
    k = np.arange(-15, 16); g = np.exp(-0.5 * (k / 5.0) ** 2); g /= g.sum()
    ind = np.stack([np.convolve(np.pad(row, 15, mode="reflect"), g, mode="valid") for row in ind])
    ind *= np.sqrt(np.mean(r0 ** 2) / np.mean(ind ** 2))
    return np.sqrt(1 - semi.ETA ** 2) * r0[:, None] + semi.ETA * ind


def reflectivity(base, tg, code, dist):
    """Deposit base[k, j] (geological time TF[k]) at the output time where tg(t, j) = TF[k] (linear split);
    samples inside a fault jump (the output step where the block code changes: missing section) are dropped;
    amplitude x0.5 within 1.5 traces of a plane."""
    nf = len(TF); R = np.zeros((nf, NX)); jump = np.diff(code, axis=0) != 0
    for j in range(NX):
        pos = np.interp(TF, tg[:, j], np.arange(nf, dtype=float))
        lo = np.clip(np.floor(pos).astype(int), 0, nf - 2)
        keep = (TF >= tg[0, j]) & (TF <= tg[-1, j]) & ~jump[lo, j]
        i = np.floor(pos).astype(int); fr = pos - i
        amp = base[:, j] * np.where(dist[np.clip(i, 0, nf - 1), j] < 1.5, 0.5, 1.0)
        for ii, ww in ((i, 1 - fr), (i + 1, fr)):
            ok = keep & (ii >= 0) & (ii < nf); np.add.at(R[:, j], ii[ok], amp[ok] * ww[ok])
    return R


def inv_q_model(cm, t):
    """1/Q of the layered model at two-way times t (derivative of lib.absorption)."""
    e = 1e-4
    return (lib.absorption(cm, t + e) - lib.absorption(cm, t - e)) / (2 * e)


def inv_q(geo, cm, c0, tg, code, j_grid):
    """True 1/Q(t, j): model Q with its step moved, x c0, x exp(layer log-Q), x block multiplier and lateral
    trend (ramped in over 1.3-1.8 s so the shallow calibration window keeps the model-like spectrum),
    then blended towards Q_GAS inside the gas lens (super-Gaussian, geological coordinates)."""
    lay = geo["layers"]; edges = np.array(lay["edges"]); vals = np.array(lay["logq"])
    k = np.clip(np.searchsorted(edges, tg, side="right") - 1, 0, len(vals) - 1)
    logq = vals[k]
    sm = np.ones(13) / 13                                               # 12 ms smoothing of the layer steps
    logq = np.apply_along_axis(lambda v: np.convolve(np.pad(v, 6, mode="edge"), sm, mode="valid"), 0, logq)
    ramp = np.clip((tg - 1.3) / 0.5, 0, 1)
    mult = np.vectorize(geo["blocks"].get)(code)
    lat = geo["trend"] * (j_grid - CENTRE) / CENTRE
    q_inv = inv_q_model(cm, tg.ravel() - lay["step_shift"]).reshape(tg.shape) / c0
    q_inv = q_inv * np.exp(-logq - ramp * (np.log(mult) + lat))
    gs = geo["gas"]
    w = np.exp(-((j_grid - gs["x"]) / gs["half_width"]) ** 4 - ((tg - gs["t"]) / gs["half_thick"]) ** 4)
    return (1 - w) * q_inv + w / gs["Q"], w


def tstar_from(qi, cm, c0):
    """Two-way t* on TF: model t* at 1.0 s / c0 above the grid, then the running integral of 1/Q."""
    t0 = lib.absorption(cm, [TF[0]])[0] / c0
    return t0 + np.vstack([np.zeros((1, qi.shape[1])), np.cumsum(0.5 * (qi[1:] + qi[:-1]) * FINE, axis=0)])


# ---------------------------------------------------------------- modelling

def synth(R, tstar, S, f, fref, N=8192):
    """Attenuated traces at TW for a t* profile per trace (tstar: (len(TF), NX)). Reflector k of trace j
    carries irfft(S AA exp(b t*_kj)); the kernel is linearly interpolated between t* grid nodes (step DTS),
    truncated to +-HALF ms. tstar=None gives the unattenuated response."""
    aa = np.clip((120 - f) / 20, 0, 1); aa = 0.5 - 0.5 * np.cos(np.pi * aa)
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / fref); b[0] = 0
    nf = len(TF); nfft = 1 << int(np.ceil(np.log2(nf + 2 * HALF + 1)))
    lagi = np.arange(-HALF, HALF + 1) % N; rows = np.round((TW - TF[0]) / FINE).astype(int) + HALF

    def kernel(ts):
        H = S * aa * (np.exp(b * ts) if ts is not None else 1); H[-1] = H[-1].real
        return np.fft.rfft(np.fft.irfft(H, N)[lagi], nfft)

    if tstar is None:
        return np.fft.irfft(np.fft.rfft(R, nfft, axis=0) * kernel(None)[:, None], nfft, axis=0)[rows]
    g = tstar / DTS; i0 = np.floor(g).astype(int); fr = g - i0
    acc = np.zeros((nfft // 2 + 1, R.shape[1]), complex)
    for gi in range(i0.min(), i0.max() + 2):
        W = R * ((i0 == gi) * (1 - fr) + (i0 == gi - 1) * fr)
        if np.any(W):
            acc += np.fft.rfft(W, nfft, axis=0) * kernel(gi * DTS)[:, None]
    return np.fft.irfft(acc, nfft, axis=0)[rows]


def fine_to_tw(A):
    return A[np.round((TW - TF[0]) / FINE).astype(int)]


def c_eff(tstar_tw, cm):
    """Effective Q factor in 1.2-3.0 s per trace: cumulative (t* ~ a_model/c) and differential (relative to
    1.2 s, what spectral ratios see); section value from the median 1/c."""
    a = lib.absorption(cm, TW); da = a - a[0]; dt = tstar_tw - tstar_tw[0]
    inv_cum = (a @ tstar_tw) / (a @ a); inv_dif = (da @ dt) / (da @ da)
    return dict(cum=1 / inv_cum, dif=1 / inv_dif, cum_section=float(1 / np.median(inv_cum)),
                dif_section=float(1 / np.median(inv_dif)))


def decay_slope(D):
    """Slope (1/Hz) of ln(coherent power 2.0-3.0 s / 1.2-1.6 s) over 10-60 Hz (descriptive realism check)."""
    f, Ps = gen.coherent_power(D, 1.2, 1.6); _, Pd = gen.coherent_power(D, 2.0, 3.0)
    m = (f >= FMAX_FIT[0]) & (f <= FMAX_FIT[1])
    return float(np.polyfit(f[m], np.log(np.maximum(Pd[m], 1e-30) / np.maximum(Ps[m], 1e-30)), 1)[0])


# ---------------------------------------------------------------- one structure

def prepare(key):
    cm = semi.curtain_cm(); lines = semi.read_lines([STRUCT[key], *NOISE[key]]); used = [int(x) for x in lines]
    D = lines[used[0]]; tau = semi.structure(semi.local_lags(D)); wv = lib.wavelets(D, TW, cm)
    noise = [semi.residual_noise(lines[L], semi.local_lags(lines[L])) for L in used[1:]]
    fc, Creal = gen.coherent_power(D, 1.2, 1.6); fq, snr_real, _ = lib.coherence_snr(D, TW, 2.0, 3.0)
    geo = draw_geology(key, np.random.default_rng(GEO_SEED[key]))
    tg, code, dist = mapping(geo, tau)
    return dict(key=key, cm=cm, D=D, tau=tau, w_in=wv["w_in"], noise=noise, lines=used, Creal=Creal,
                band=(fc >= 5) & (fc <= 80), snr_real=snr_real, bq=(fq >= 10) & (fq <= 80), geo=geo,
                tg=tg, code=code, dist=dist)


def calibrate_level(ctx, x, N0):
    bq = ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + L * N0, TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in gen.LEVEL_GRID]
    i = int(np.argmin(err)); return float(gen.LEVEL_GRID[i]), err[i]


def build(ctx, c0, r):
    """Z and Y datasets for one structure, Q level c0 and realisation r (reflectivity seed, noise line r)."""
    cm, key = ctx["cm"], ctx["key"]; f = np.fft.rfftfreq(8192, FINE); fref = float(cm["fref"])
    base = base_series(np.random.default_rng(RC_SEED[(key, r)]))
    R = reflectivity(base, ctx["tg"], ctx["code"], ctx["dist"])
    jg = np.broadcast_to(np.arange(NX), ctx["tg"].shape)
    qi, gas_w = inv_q(ctx["geo"], cm, c0, ctx["tg"], ctx["code"], jg)
    ts_z = tstar_from(qi, cm, c0); ce = c_eff(fine_to_tw(ts_z), cm)
    ts_y = np.repeat(semi.tstar_true(cm, ce["cum_section"])[:, None], NX, axis=1)
    S0 = semi.source_spectrum(ctx["w_in"], f)
    a14 = lib.absorption(cm, [1.4])[0]; out = {}
    k_scale = level = None; N0 = ctx["noise"][r]
    for var, ts in (("Z", ts_z), ("Y", ts_y)):
        S = S0 * np.exp(np.pi * f * (np.median(np.interp(1.4, TF, ts[:, CENTRE])) - a14))
        x = synth(R, ts, S, f, fref)
        if var == "Z":                                                   # Y reuses Z's scale and noise level
            b = ctx["band"]; k_scale = float(np.sqrt(ctx["Creal"][b].sum() / semi.band_power(x, 1.2, 1.6)[b].sum()))
            level, cal_err = calibrate_level(ctx, x * k_scale, N0)
        x = x * k_scale; noise = level * N0; f_true = semi.true_band(x, noise)
        src = synth(R, None, S * k_scale, f, fref)
        ideal = synth(R, None, S * k_scale * semi.cut(f, f_true), f, fref)
        D = x + noise
        meta = dict(id=f"{var}{key}_c{c0}_r{r}", variant=var, structure=key, inline=ctx["lines"][0],
                    noise_line=ctx["lines"][1 + r], c0=c0, r=r, rc_seed=RC_SEED[(key, r)], scale=k_scale,
                    level=level, f_true=float(f_true), decay_slope=decay_slope(D),
                    c_eff_cum=float(c_eff(fine_to_tw(ts), cm)["cum_section"]),
                    c_eff_dif=float(c_eff(fine_to_tw(ts), cm)["dif_section"]))
        if var == "Z":
            meta.update(snr_cal_err=cal_err, c_eff_cum_trace=[float(v) for v in ce["cum"]],
                        c_eff_dif_trace=[float(v) for v in ce["dif"]])
        out[var] = dict(meta=meta, arrays=dict(
            D=D, signal=x, noise=noise, ideal=ideal, src=src, R=R,
            tstar=fine_to_tw(ts), invQ=fine_to_tw(qi) if var == "Z" else fine_to_tw(np.gradient(ts, FINE, axis=0)),
            block=fine_to_tw(ctx["code"]), fault_dist=fine_to_tw(ctx["dist"]),
            gas=fine_to_tw(gas_w) if var == "Z" else np.zeros((len(TW), NX))))
    return out


def load(did):
    """Arrays (float32, TW x NX) and metadata of a built dataset, e.g. load("ZA_c1.05_r0")."""
    root = semi.ROOT.parent / "semisynth_fault_q_20261030"
    meta = json.loads((root / "analysis" / "datasets.json").read_text(encoding="utf-8"))[did]
    return dict(np.load(root / "cache" / "data" / f"{did}.npz")), meta


def jsonable(geo):
    g = dict(geo); g["blocks"] = {str(k): v for k, v in geo["blocks"].items()}; return g


if __name__ == "__main__":
    print(json.dumps(jsonable(draw_geology("A", np.random.default_rng(GEO_SEED["A"]))), indent=1))
