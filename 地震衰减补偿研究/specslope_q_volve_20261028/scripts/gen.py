"""Semi-synthetic datasets of experiment 10 with the experiment-8 recipe (semi.py, c395b78): F-11 T2 log
reflectivity laid along a real inline's structure, smooth true Q = c x model Q, source from the real shallow
spectrum, real noise from other inlines (dip-aligned lateral second difference), level calibrated so that
the 2.0-3.0 s coherent SNR matches the real line in 10-80 Hz.

prepare(inline, noise_lines) -> context; signal(ctx, c, seed) -> dict(scaled signal, calibrated noise level, reflectivity, scale, source).
With segs (wells.py segments) the reflectivity is the colour stress version: whole real well segments (wells.py) in random
order instead of the 40 ms block bootstrap of the F-11 T2 log, so that windows carry the real well-to-well
and depth-to-depth reflectivity colour.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq  # noqa: E402
from sq import np, lib, semi, TW  # noqa: E402

LEVEL_GRID = np.exp(np.linspace(np.log(0.5), np.log(4), 43))


def coherent_power(D, a, b):
    f, snr, P = lib.coherence_snr(D, TW, a, b)
    return f, P * snr / (1 + snr)


def prepare(inline, noise_lines):
    cm = semi.curtain_cm(); lines = semi.read_lines([inline] + list(noise_lines)); used = list(lines)
    D = lines[used[0]]; lags = semi.local_lags(D); tau = semi.structure(lags); wv = lib.wavelets(D, TW, cm)
    noise = [semi.residual_noise(lines[L], semi.local_lags(lines[L])) for L in used[1:]]
    fc, Creal = coherent_power(D, 1.2, 1.6); fq, snr_real, _ = lib.coherence_snr(D, TW, 2.0, 3.0)
    return dict(cm=cm, D=D, tau=tau, w_in=wv["w_in"], noise=noise, lines=[int(x) for x in used],
                Creal=Creal, band=(fc >= 5) & (fc <= 80), snr_real=snr_real, bq=(fq >= 10) & (fq <= 80))


def well_reflectivity(tau, rng, segs):
    """As semi.reflectivity, but every series is a concatenation of whole well segments (random order; the first
    one starts at a random offset); the centre series is one such concatenation."""
    nf, nx = len(semi.TF), tau.shape[1]

    def series(n):
        out, first = [], True
        while sum(len(x) for x in out) < n:
            s = segs[rng.integers(len(segs))]
            out.append(s[rng.integers(len(s)):] if first else s); first = False
        return np.concatenate(out)[:n]

    r0 = series(nf); ind = np.stack([series(nf) for _ in range(nx)], axis=1)
    k = np.arange(-15, 16); g = np.exp(-0.5 * (k / 5.0) ** 2); g /= g.sum()
    ind = np.stack([np.convolve(np.pad(row, 15, mode="reflect"), g, mode="valid") for row in ind])
    ind *= np.sqrt(np.mean(r0 ** 2) / np.mean(ind ** 2))
    base = np.sqrt(1 - semi.ETA ** 2) * r0[:, None] + semi.ETA * ind
    R = np.zeros((nf, nx))
    for j in range(nx):
        pos = (semi.TF + np.interp(semi.TF, TW, tau[:, j]) - semi.TF[0]) / semi.FINE
        i = np.floor(pos).astype(int); fr = pos - i
        for ii, ww in ((i, 1 - fr), (i + 1, fr)):
            ok = (ii >= 0) & (ii < nf); np.add.at(R[:, j], ii[ok], base[ok, j] * ww[ok])
    return R


def signal(ctx, c, seed, segs=None):
    """Signal for true factor c and reflectivity seed (well segments if segs is given); noise level
    calibrated against noise[0]. Returns dict(x, level, R, scale, S)."""
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE); rng = np.random.default_rng(seed)
    R = semi.reflectivity(ctx["tau"], rng)[0] if segs is None else well_reflectivity(ctx["tau"], rng, segs)
    S = semi.source_spectrum(ctx["w_in"], f, 0.0, cm, c)
    x = semi.model(R, semi.tstar_true(cm, c), S, f, float(cm["fref"]))
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / semi.band_power(x, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + L * N0, TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in LEVEL_GRID]
    return dict(x=x, level=float(LEVEL_GRID[int(np.argmin(err))]), R=R, scale=k, S=S, f=f)
