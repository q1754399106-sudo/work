"""Shared imports of uqw_snr_20261114 (development data only): the UQ-W improvement code of uqw_improve_20261113 unchanged."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "uqw_improve_20261113" / "scripts"))
from probe2 import np, semi, score, sc11, TW, DT, snr_rows, CW, load, lib, t2, lib12, SEL, tvfilt, CURVE, MS  # noqa: E402,F401
import uqwh as H  # noqa: E402,F401
from stage0 import DEV  # noqa: E402,F401

fd = CW.fd
WIN, HOP, NF = 0.4, 0.1, 512                      # the windows of uqw.snr_rows / lib.coherence_snr
F5 = np.fft.rfftfreq(NF, DT); F10 = np.fft.rfftfreq(1024, DT)
BANDS = ((10, 25), (25, 40), (40, 55), (55, 70))


def noise_of(arr, kind, x, seed):
    """Noise field of a development set: 'res' the real residual noise (lateral second difference), 'gau' spatially white
    noise with the same time-variant spectrum and rms (semi.gaussian_noise); x the noise multiplier."""
    return x * (arr["noise"] if kind == "res" else semi.gaussian_noise(arr["noise"], np.random.default_rng(seed)))


def spectra(X, tw, a, b):
    w = (tw >= a - 1e-9) & (tw <= b + 1e-9)
    return np.fft.rfft(X[w] * lib.taper(w.sum())[:, None], NF, axis=0)


def cross(Fx, k, sm=5):
    """Mean real cross-power of traces k apart (k = 0: mean power), 5-bin smoothing as lib.coherence_snr."""
    if k == 0:
        return lib.smooth(np.mean(np.abs(Fx) ** 2, axis=1), sm)
    return lib.smooth(np.mean(np.real(Fx[:, :-k] * np.conj(Fx[:, k:])), axis=1), sm)


def starts(tw):
    return np.arange(tw[0], tw[-1] - WIN + 1e-9, HOP)
