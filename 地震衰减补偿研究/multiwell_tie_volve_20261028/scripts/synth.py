"""Semi-synthetic wells for the calibration of WR (experiment-8/10 recipe with a real well log at the centre column).

Earth reflectivity of the centre ("well") trace: the well's Gardner RC at its own times, warped by a smooth random
time error (log-to-seismic time mismatch, sd JIT s, correlation length 0.3 s) and mixed with an independent
series, earth = sqrt(1 - eps^2) * log + eps * other (log errors, upscaling, 3-D effects); outside the log span a
circular 40 ms block bootstrap of the same log. Neighbouring traces as semi.reflectivity (eta = 0.3, laterally
smoothed independent part), displaced along the real inline's structure. Attenuation with true factor c, source,
scaling and real-noise level calibration exactly as gen.signal (experiment 10). The estimator only sees the log RC.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq  # noqa: E402
from wq import np, semi, gen  # noqa: E402

JIT = 0.008
EPS = 0.3
BLOCK = 40


def warp(rng, n, sd=JIT, corr=0.3):
    """Smooth random time error (s) on a 1 ms grid: Gaussian-filtered white noise scaled to sd."""
    k = np.arange(-int(2 * corr / wq.FINE), int(2 * corr / wq.FINE) + 1); g = np.exp(-0.5 * (k * wq.FINE / (corr / 2)) ** 2)
    e = np.convolve(rng.standard_normal(n + len(k)), g, mode="same")[len(k) // 2:len(k) // 2 + n]
    return e / e.std() * sd


def centre_series(w, rng, eps=EPS, jit=JIT):
    """(log series on semi.TF, earth series on semi.TF)."""
    nf = len(semi.TF); rc = w["rc"]

    def boot(n):
        s = rng.integers(0, len(rc), int(np.ceil(n / BLOCK)))
        return np.concatenate([np.take(rc, np.arange(k, k + BLOCK), mode="wrap") for k in s])[:n]

    i0 = int(round((w["t"][0] - semi.TF[0]) / wq.FINE)); i1 = min(i0 + len(rc), nf)
    log = boot(nf); log[i0:i1] = rc[:i1 - i0]
    shifted = np.interp(semi.TF + warp(rng, nf, jit), semi.TF, log)
    return log, np.sqrt(1 - eps ** 2) * shifted + eps * boot(nf)


def reflectivity(ctx, w, rng, eps=EPS, jit=JIT):
    tau = ctx["tau"]; nf, nx = len(semi.TF), tau.shape[1]; rc = w["rc"]
    log, r0 = centre_series(w, rng, eps, jit)

    def boot(n):
        s = rng.integers(0, len(rc), int(np.ceil(n / BLOCK)))
        return np.concatenate([np.take(rc, np.arange(k, k + BLOCK), mode="wrap") for k in s])[:n]

    ind = np.stack([boot(nf) for _ in range(nx)], axis=1)
    k = np.arange(-15, 16); g = np.exp(-0.5 * (k / 5.0) ** 2); g /= g.sum()
    ind = np.stack([np.convolve(np.pad(row, 15, mode="reflect"), g, mode="valid") for row in ind])
    ind *= np.sqrt(np.mean(r0 ** 2) / np.mean(ind ** 2))
    base = np.sqrt(1 - semi.ETA ** 2) * r0[:, None] + semi.ETA * ind
    base[:, semi.CENTRE] = r0                       # the well trace carries the earth series itself
    R = np.zeros((nf, nx))
    for j in range(nx):
        pos = (semi.TF + np.interp(semi.TF, semi.TW, tau[:, j]) - semi.TF[0]) / wq.FINE
        i = np.floor(pos).astype(int); fr = pos - i
        for ii, ww in ((i, 1 - fr), (i + 1, fr)):
            ok = (ii >= 0) & (ii < nf); np.add.at(R[:, j], ii[ok], base[ok, j] * ww[ok])
    return R


def signal(ctx, w, c, seed, eps=EPS, jit=JIT):
    """gen.signal with the well reflectivity: dict(x, level)."""
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE); rng = np.random.default_rng(seed)
    R = reflectivity(ctx, w, rng, eps, jit)
    S = semi.source_spectrum(ctx["w_in"], f, 0.0, cm, c)
    x = semi.model(R, semi.tstar_true(cm, c), S, f, float(cm["fref"]))
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / semi.band_power(x, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(wq.lib.coherence_snr(x + L * N0, wq.TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in gen.LEVEL_GRID]
    return dict(x=x, level=float(gen.LEVEL_GRID[int(np.argmin(err))]))


def sections(D, wins, half=None):
    j, h = semi.CENTRE, wq.HALF if half is None else half
    return [D[:, max(j - h, 0):j + h + 1] for _ in wins]
