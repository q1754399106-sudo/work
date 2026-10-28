"""Experiment 8 runs (PLAN_ZH.md).

usage: run.py build              build signals, noises and the realism check (no method is run)
       run.py solve K N          worker K of N: run every method on its datasets (resumable)
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import semi  # noqa: E402
from semi import np, lib, bm, TW, TF, C_TRUE  # noqa: E402

ROOT = semi.ROOT
CACHE, AN, SC = ROOT / "cache", ROOT / "analysis", ROOT / "results" / "scores"
BUILD = CACHE / "build.npz"
PHASE = 30.0


def datasets():
    """(id, signal key, noise key, noise level): primary first, then the three secondary groups. Levels are
    the calibrated ones of analysis/build.json, except group R1 (residual noise at its own level)."""
    lv = json.loads((AN / "build.json").read_text(encoding="utf-8"))["level"]
    out = [(f"P_c{c}_n{i}", f"c{c}", f"R{i}", lv[f"c{c}"]) for c in C_TRUE for i in range(10)]
    out += [(f"G_c{c}_n{i}", f"c{c}", f"G{i}", lv[f"c{c}"]) for c in C_TRUE for i in range(5)]
    out += [(f"R1_c{c}_n{i}", f"c{c}", f"R{i}", 1.0) for c in C_TRUE for i in range(5)]
    out += [(f"PH_c1.0_n{i}", "ph30", f"R{i}", lv["ph30"]) for i in range(5)]
    return out


def haar_sigma(x):
    d = (x[1::2][: len(x) // 2] - x[0::2][: len(x) // 2]) / np.sqrt(2)
    return float(np.median(np.abs(d)) / 0.6745)


def coherent_power(D, a, b):
    f, snr, P = lib.coherence_snr(D, TW, a, b)
    return f, P * snr / (1 + snr)


def build():
    t0 = time.monotonic(); cm = semi.curtain_cm()
    lines = semi.read_lines([semi.INLINE] + semi.NOISE_LINES)
    used = list(lines); D = lines[used[0]]; assert used[0] == semi.INLINE
    lags = semi.local_lags(D); tau = semi.structure(lags)
    wv = lib.wavelets(D, TW, cm)
    rng = np.random.default_rng(semi.SEED)
    R, r0 = semi.reflectivity(tau, rng)
    f = np.fft.rfftfreq(8192, semi.FINE); fref = float(cm["fref"])
    fc, Creal = coherent_power(D, 1.2, 1.6); band = (fc >= 5) & (fc <= 80)
    noise = {}
    for i, L in enumerate(used[1:]):
        noise[f"R{i}"] = semi.residual_noise(lines[L], semi.local_lags(lines[L]))
    E0 = semi.residual_noise(D, lags)
    for s in range(5):
        noise[f"G{s}"] = semi.gaussian_noise(E0, np.random.default_rng(semi.SEED + 1 + s))
    fq, snr_real, _ = lib.coherence_snr(D, TW, 2.0, 3.0); bq = (fq >= 10) & (fq <= 80)
    grid = np.exp(np.linspace(np.log(0.5), np.log(4), 43))
    sig, scale, level, misfit, src = {}, {}, {}, {}, {}
    for key, c, ph in [(f"c{c}", c, 0.0) for c in C_TRUE] + [("ph30", 1.0, PHASE)]:
        S = semi.source_spectrum(wv["w_in"], f, ph, cm, c)
        x = semi.model(R, semi.tstar_true(cm, c), S, f, fref)
        k = float(np.sqrt(Creal[band].sum() / semi.band_power(x, 1.2, 1.6)[band].sum()))
        sig[key] = x * k; scale[key] = k
        if ph == 0:
            src[key] = semi.model(R, None, semi.source_spectrum(wv["w_in"], f, 0.0, cm, c), f, fref) * k
        err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(sig[key] + L * noise["R0"], TW, 2.0, 3.0)[1][bq], 1e-3))
                                      - np.log(np.maximum(snr_real[bq], 1e-3))) ** 2))) for L in grid]
        level[key] = float(grid[int(np.argmin(err))]); misfit[key] = float(min(err))
    np.savez_compressed(BUILD, R=R.astype(np.float32), r0=r0, tau=tau, w_in=wv["w_in"], E0=E0, D_real=D,
                        **{f"src_{k}": v for k, v in src.items()},
                        **{f"sig_{k}": v for k, v in sig.items()}, **{f"noise_{k}": v for k, v in noise.items()},
                        **{f"scale_{k}": v for k, v in scale.items()})
    # Realism check: blind wavelet/band estimates and the target/shallow spectral decay, real vs synthetic.
    def decay(X):
        fa, _, Ps = lib.coherence_snr(X, TW, 1.2, 1.6); _, _, Pd = lib.coherence_snr(X, TW, 2.0, 3.0)
        b = (fa >= 10) & (fa <= 60)
        return fa[b], 0.5 * np.log(Pd[b] / Ps[b])
    fr, dr = decay(D); real = dict(f_peak=wv["f_peak"], fmax_shallow=wv["fmax_shallow"], fmax_target=wv["fmax_target"])
    slope = lambda fx, dx: float(np.polyfit(fx, dx, 1)[0])  # noqa: E731
    real["decay_slope_per_Hz"] = slope(fr, dr)
    rows = {}
    for key in sig:
        for nk in ("R0", "G0"):
            X = sig[key] + level[key] * noise[nk]; w = lib.wavelets(X, TW, cm); fx, dx = decay(X)
            rows[f"{key}+{nk}"] = dict(f_peak=w["f_peak"], fmax_shallow=w["fmax_shallow"], fmax_target=w["fmax_target"],
                                      decay_slope_per_Hz=slope(fx, dx), rms_log_decay_vs_real=float(np.sqrt(np.mean((dx - dr) ** 2))),
                                      f_true=semi.true_band(sig[key], level[key] * noise[nk]))
    _, Cn = coherent_power(D, 2.0, 3.0); _, _, Pn = lib.coherence_snr(D, TW, 2.0, 3.0)
    Pe = semi.band_power(E0, 2.0, 3.0); b = (fc >= 10) & (fc <= 80)
    nrat = float(np.sqrt(Pe[b].sum() / (Pn - Cn)[b].sum()))
    bm.writejson(AN / "build.json", dict(lines=[int(x) for x in used], noise_lines=[int(x) for x in used[1:]],
                                         scale=scale, level=level, snr_log_misfit=misfit, real=real, synthetic=rows,
                                         residual_vs_coherence_noise_rms_10_80Hz=nrat,
                                         lag_abs_median_samples=float(np.median(np.abs(lags))),
                                         tau_range_ms=[float(tau.min() * 1000), float(tau.max() * 1000)],
                                         seconds=time.monotonic() - t0))
    print(json.dumps(dict(real=real, level=level, misfit=misfit, synthetic=rows, noise_ratio=nrat), indent=1), flush=True)


def methods(D, cm, c):
    """Every method of PLAN_ZH.md on one dataset; outputs and diagnostics."""
    live = np.abs(D).sum(axis=0) > 0; K, _, _ = bm.predictor(D, live)
    sigma = np.array([haar_sigma(D[:, j]) for j in range(D.shape[1])])
    Q = np.asarray(cm["Q"], float); wv = lib.wavelets(D, TW, cm); fc = min(wv["fmax_target"], wv["fmax_shallow"])
    out, info = {"Input": D}, {}
    for g in (10, 20, 30):
        out[f"IQ{g}"] = semi.inverse_q(D, cm, Q, g)
    out["IQ20a"] = semi.lowpass(out["IQ20"], fc)
    B, _, P, scale = bm.dictionary(cm)
    out["SC"], info["SC"] = bm.solve(B, P, D, sigma, scale, K, *lib.SEL_SC, live)
    out["SCa"], info["SCa"], _ = lib.solve(cm, Q, wv, D, sigma, K, live, *lib.SEL_SC)
    factor, qi = lib.qscan(cm, Q, wv, D, sigma, K, live); info["qscan"] = dict(factor=factor, **qi)
    out["QSCa"], info["QSCa"], _ = lib.solve(cm, Q * factor, wv, D, sigma, K, live, *lib.SEL_SC)
    out["QIQa"] = semi.lowpass(semi.inverse_q(D, cm, Q * factor, 20), fc)
    out["SCa_oracle"], info["SCa_oracle"], _ = lib.solve(cm, Q * c, wv, D, sigma, K, live, *lib.SEL_SC)
    info["wavelet"] = dict(f_peak=wv["f_peak"], fmax_shallow=wv["fmax_shallow"], fmax_target=wv["fmax_target"], fc=fc)
    return out, info


def solve(k, n):
    cm = semi.curtain_cm(); B = np.load(BUILD); f = np.fft.rfftfreq(8192, semi.FINE)
    R = B["R"].astype(float)
    tie_ref = semi.ricker_ref(R[:, semi.CENTRE], float(cm["fpeak"]))
    (CACHE / "out").mkdir(exist_ok=True); SC.mkdir(parents=True, exist_ok=True)
    for j, (did, skey, nkey, level) in enumerate(datasets()):
        if j % n != k or (SC / f"{did}.json").exists():
            continue
        t0 = time.monotonic(); c = 1.0 if skey == "ph30" else float(skey[1:]); ck = f"c{c}"
        sig, noise = B[f"sig_{skey}"], level * B[f"noise_{nkey}"]
        D = sig + noise; ftrue = semi.true_band(sig, noise); src = B[f"src_{ck}"]
        S = semi.source_spectrum(B["w_in"], f, 0.0, cm, c) * float(B[f"scale_{ck}"])
        ideal = semi.model(R, None, S * semi.cut(f, ftrue), f, float(cm["fref"]))
        out, info = methods(D, cm, c)
        sc = {m: semi.scores(A, ideal, src, tie_ref) for m, A in out.items()}
        bm.writejson(SC / f"{did}.json", dict(id=did, c_true=c, signal=skey, noise=nkey, level=level, f_true=ftrue,
                                              scores=sc, info=info, seconds=time.monotonic() - t0))
        np.savez_compressed(CACHE / "out" / f"{did}.npz", ideal=ideal.astype(np.float32),
                            **{m: A.astype(np.float32) for m, A in out.items()})
        print(did, f"{time.monotonic() - t0:.0f}s qscan {info['qscan']['factor']:.3f} ftrue {ftrue:.1f}",
              {m: round(sc[m]['M1'], 3) for m in ('Input', 'IQ20a', 'SCa', 'QSCa')}, flush=True)


if __name__ == "__main__":
    AN.mkdir(exist_ok=True); CACHE.mkdir(exist_ok=True)
    if sys.argv[1] == "build":
        build()
    else:
        solve(int(sys.argv[2]), int(sys.argv[3]))
