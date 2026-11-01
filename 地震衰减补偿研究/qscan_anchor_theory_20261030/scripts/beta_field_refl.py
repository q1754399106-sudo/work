"""Exploratory (BETA_ZH.md step 3): semi-synthetic data built on the field's own reflectivity estimate.

usage: beta_field_refl.py    worker; jobs claimed with results13/<job>.claim, resumable.
R3  reflectivity estimated from a real inline (1.0-3.2 s): stabilised inverse Q at its own SR factor (70 dB, as
    lib12.iq), then division by the input-wavelet spectrum |w_in| inside a band mask (6 Hz to the coherent-SNR
    limit of the line, cosine tapers), upsampled to 1 ms. Lines 10193, 10241, 10305 (SR lines, not used as
    structure or noise lines here).
R0b control: the experiment-14 reflectivity (gen.signal's R, structure 10227, seeds 9995-9997) passed through the
    same band mask (of the line with the same index).
Both are forward modelled exactly as gen.signal (source of the structure line 10227, true Q = c x curtain Q,
shallow coherent power matched, noise line 10149 with the calibrated level x 1), c 0.8 / 1.0 / 1.25, and scanned
with fixed anchors 0.7 / 1.0 / 1.4 (lambda 1/2, qprofile.scan).
-> results13/<variant>_<k>_c<c>__a<c_a>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from qprofile import scan

sys.path.insert(0, str(ROOT.parent / "field_baselines_volve_20261029" / "scripts"))
import lib13 as L13  # noqa: E402

OUT = ROOT / "results13"; bm.MAXITER = 20000
LINES = (10193, 10241, 10305); SEEDS = (9995, 9996, 9997); CS = (0.8, 1.0, 1.25); ANCH = (0.7, 1.0, 1.4)
TWIDE = np.round(np.arange(int(round(1.0 / semi.DT)), int(round(3.2 / semi.DT)) + 1) * semi.DT, 6)
SRL = json.loads((L13.IN / "specslope_q_volve_20261028" / "results" / "real" / "lines.json").read_text())


def read_wide(L):
    h, x, y, il, xl = L13.cube(); w = semi.win(np.arange(1126) * semi.DT, 1.0, 3.2)
    D = np.zeros((int(w.sum()), (L13.XL1 - L13.XL0) // L13.XSTEP + 1))
    for j in np.where((il == L) & (xl >= L13.XL0) & (xl <= L13.XL1) & ((xl - L13.XL0) % L13.XSTEP == 0))[0]:
        D[:, (xl[j] - L13.XL0) // L13.XSTEP] = semi.ibm(h[j, 240:].copy().view(">u4"))[w]
    return D


def inverse_q(D, cm, Q, c, g_db=70, nfft=2048):
    """lib12.iq on the wide time axis TWIDE (amplitude and dispersion at factor c)."""
    n = len(TWIDE); f = np.fft.rfftfreq(nfft, semi.DT); lag = np.arange(n)
    t = lib.absorption(dict(cm, Q=Q * c), TWIDE); s = semi.gain_sigma2(g_db)
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); M = np.empty((n, n))
    for i in range(n):
        b = np.exp(-np.pi * f * t[i]); H = np.exp(-2j * f * t[i] * lf) * (b + s) / (b ** 2 + s)
        H[0] = H[0].real; H[-1] = H[-1].real; M[i] = np.fft.irfft(H, nfft)[(i - lag) % nfft]
    return M @ D


def band_mask(f, lo, hi):
    m = np.ones_like(f); r = (f > lo - 3) & (f < lo); m[r] = 0.5 - 0.5 * np.cos(np.pi * (f[r] - (lo - 3)) / 3)
    m[f <= lo - 3] = 0; r = (f > hi) & (f < hi + 8); m[r] = 0.5 + 0.5 * np.cos(np.pi * (f[r] - hi) / 8); m[f >= hi + 8] = 0
    return m


def field_reflectivity(L, cm, Q, w_in):
    D = read_wide(L); live = np.abs(D).sum(axis=0) > 0; c = float(SRL[str(L)]["c"])
    X = inverse_q(D, cm, Q, c)
    _, snr, _ = lib.coherence_snr(D[(TWIDE >= 1.2) & (TWIDE <= 3.0 + 1e-9)], semi.TW, 2.0, 3.0)
    f4, _, P = lib.coherence_snr(D[(TWIDE >= 1.2) & (TWIDE <= 3.0 + 1e-9)], semi.TW, 1.2, 1.6)
    fpk = float(f4[np.argmax(np.sqrt(P) * (f4 >= 5))]); hi = float(min(lib.snr_band(f4, snr, fpk), 80.0))
    N = 2048; f = np.fft.rfftfreq(N, semi.DT); W = np.abs(np.fft.rfft(np.roll(np.pad(w_in, (0, N - len(w_in))), -100), N))
    m = band_mask(f, 6.0, hi); inv = np.where(W > 1e-3 * W.max(), m / np.maximum(W, 1e-12), 0) * W.max()
    R4 = np.fft.irfft(np.fft.rfft(X, N, axis=0) * inv[:, None], N, axis=0)[:len(TWIDE)]
    # to 1 ms on semi.TF (band-limited, FFT interpolation)
    n1 = len(semi.TF); Rf = np.fft.rfft(R4 * live[None, :], axis=0)
    R1 = np.fft.irfft(Rf, len(TWIDE) * 4, axis=0)[:n1] * 4
    return R1, hi, live


def mask_1ms(R, hi):
    N = 8192; f = np.fft.rfftfreq(N, semi.FINE); m = band_mask(f, 6.0, hi)
    return np.fft.irfft(np.fft.rfft(R, N, axis=0) * m[:, None], N, axis=0)[:R.shape[0]]


def signal_R(ctx, c, R):
    """gen.signal with a given 1 ms reflectivity R (len(semi.TF), nx)."""
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE)
    S = semi.source_spectrum(ctx["w_in"], f, 0.0, cm, c); x = semi.model(R, semi.tstar_true(cm, c), S, f, float(cm["fref"]))
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / gen.semi.band_power(x, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + Lv * N0, semi.TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for Lv in gen.LEVEL_GRID]
    return dict(x=x, level=float(gen.LEVEL_GRID[int(np.argmin(err))]))


def jobs():
    return [(v, k, c, a) for k in range(3) for v in ("R3", "R0b") for c in CS for a in ANCH]


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(10227, [10149, 10273, 10355]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    cacheR = {}
    for v, k, c, a in jobs():
        job = f"{v}_{k}_c{c}__a{a}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic()
        if ("F", k) not in cacheR:
            cacheR[("F", k)] = field_reflectivity(LINES[k], cm, Q, ctx["w_in"])
        R1, hi, live = cacheR[("F", k)]
        if v == "R3":
            R = R1[:, :ctx["tau"].shape[1]]
        else:
            if ("B", k, c) not in cacheR:
                cacheR[("B", k, c)] = mask_1ms(gen.signal(ctx, c, SEEDS[k])["R"], hi)
            R = cacheR[("B", k, c)]
        g = signal_R(ctx, c, R); D = g["x"] + g["level"] * ctx["noise"][0]; pp = sc11.prepare_sc(D, cm)
        wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * a)); res = scan(D, cm, Q, pp, lambda cc: wv)
        rec = dict(variant=v, k=k, line=LINES[k], seed=SEEDS[k], band_hi=hi, c_true=c, u=1 / c, c_anchor=a, v_anchor=1 / a,
                   level=g["level"], seconds=time.monotonic() - t0, **res)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", f"hi {hi:.0f}", f"{rec['seconds']:.0f}s", flush=True)
