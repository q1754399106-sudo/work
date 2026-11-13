"""UQ-W followed by half-strength spectral shaping (candidate "UQW-H", NOTE_ZH.md section 6).

uqwh(D, cm, Q, c, tw, anchor, w_out, inv) = half(uqw.uqw(D, cm, Q, c, 0.447, tw, anchor, 0, inv_snr=inv), w_out)
half: in each Hann window (0.5 s, hop 0.1 s) the filter of the VDB shaping, G0 = |W_out| A / (A^2 + eps max A^2) with A the
section-mean amplitude spectrum (11-bin smoothing), is applied as G0^LAM inside LO-HI Hz and as G0 outside. LAM = 0.5.
Same statements as dev3.shape(kind="lam"); the time axis enters only through the number of samples (4 ms sampling).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "qunknown_comp_20261109" / "scripts"))
from uqw import np, uqw, DT  # noqa: E402
from stage0 import lib, SEL  # noqa: E402

LAM, LO, HI, NF, WIN, HOP, NSM, SIG, GAMMA = 0.5, 8.0, 72.0, 512, SEL["B3"]["win"], 0.1, 11, 0.447, 0.0


def half(X, w_out, lam=LAM):
    n = X.shape[0]; L = int(round(WIN / DT)); step = int(round(HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(NF, DT)
    live = np.abs(X).sum(axis=0) > 0; Wo = np.abs(np.fft.rfft(w_out, NF)); Wo = Wo / Wo.max(); band = (f >= LO) & (f <= HI)
    out = np.zeros_like(X); ws = np.zeros(n)
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; Fx = np.fft.rfft(X[i0:i1] * hw[:, None], NF, axis=0)
        A = np.sqrt(lib.smooth(np.mean(np.abs(Fx[:, live]) ** 2, axis=1), NSM)); A = A / A[band].max()
        G0 = Wo * A / (A ** 2 + SEL["B3"]["eps"] * np.max(A ** 2)); G = np.where(band, G0 ** lam, G0)
        out[i0:i1] += np.fft.irfft(Fx * G[:, None], NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]


def uqwh(D, cm, Q, c, tw, anchor, w_out, inv):
    return half(uqw(D, cm, Q, c, SIG, tw, anchor, GAMMA, inv_snr=inv), w_out)
