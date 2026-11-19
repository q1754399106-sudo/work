"""Frozen time-invariant shaping of PLAN_CONFIRM_SHAPING_ZH.md section 2.1 (used by UQW_TI and PH_TI).

target(U)        T(f): section-mean amplitude spectrum of the UQ-W output U in the anchor window 1.2-1.6 s (live traces, Hann,
                 nfft 1024, power smoothed over 21 bins ~ 5 Hz, square root), normalised to its maximum in 8-72 Hz.
ti_shape(X, T)   A(f): the same spectrum of X over 1.2-3.0 s; G = T / max(A, 1e-3) x taper (dev3.taper: 1 in 8-72 Hz, cosine
                 to 0 over 4 Hz below and 12 Hz above); output irfft(rfft(X, 1024) G)[:n]. Full strength, one filter per section.
Inputs are output sections and the time axis only; nothing else is read.
"""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BASE / "uqw_improve_20261113" / "scripts"))
from dev3 import taper, np, lib  # noqa: E402

NFFT, NSM, ANCHOR, EST, LO, HI, FLOOR = 1024, 21, (1.2, 1.6), (1.2, 3.0), 8.0, 72.0, 1e-3


def _amp(X, tw, win, f):
    k = (tw >= win[0]) & (tw <= win[1]); live = np.abs(X).sum(axis=0) > 0
    P = np.mean(np.abs(np.fft.rfft(X[k][:, live] * np.hanning(k.sum())[:, None], NFFT, axis=0)) ** 2, axis=1)
    A = np.sqrt(lib.smooth(P, NSM)); band = (f >= LO) & (f <= HI)
    return A / A[band].max()


def target(U, tw, dt):
    return _amp(U, tw, ANCHOR, np.fft.rfftfreq(NFFT, dt))


def ti_shape(X, T, tw, dt):
    n = X.shape[0]; f = np.fft.rfftfreq(NFFT, dt); G = T / np.maximum(_amp(X, tw, EST, f), FLOOR) * taper(f)
    return np.fft.irfft(np.fft.rfft(X, NFFT, axis=0) * G[:, None], NFFT, axis=0)[:n]
