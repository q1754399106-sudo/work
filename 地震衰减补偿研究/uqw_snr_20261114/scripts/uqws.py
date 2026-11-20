"""UQW-S: UQ-W followed by half-strength shaping, restricted to the band where the target wavelet is measured.

Target (no free parameter added):
  w_t = w_in of sca_ceiling libc.wavelets(robust=True): the shallow (1.2-1.6 s) coherent wavelet, de-attenuated to the anchor,
        cut at the gap-tolerant shallow SNR edge fs_r; no deep-window cut (that cut collapses under noise, dev6/dev9).
  T(t, f) = |W_t(f)| sqrt(w(t, f)), w = snr / (1 + snr) (window mean of UQ-W's coherence snr): the expected output amplitude
        of a Wiener filter.
Filter per Hann window (0.5 s, hop 0.1 s, as uqwh.half): G0 = T A / (A^2 + eps max A^2), A the output amplitude spectrum,
  G = G0^LAM for LO <= f <= fs_r - 5 Hz, blended (cosine over 5 Hz) to 1 at fs_r (above it UQ-W is left alone: no target there),
  and G0 above HI = 72 Hz (removes the out-of-band gain of UQ-W found on the field data). LAM = 0.5 as uqwh.py.
fs_r = inf (no cap) when the target is known (FD: Ricker).
"""
import sys
from env import np, lib, DT, H, F10, SEL, BASE, TW
sys.path.insert(0, str(BASE / "fault_q_confirm_20261030" / "scripts")); sys.path.insert(0, str(BASE / "sca_ceiling_20261030" / "scripts"))
import libc  # noqa: E402

LAM, TAPER = 0.5, 5.0


def target(D, cm):
    wr = libc.wavelets(D, cm, True); return wr["w_in"], float(wr["fmax_shallow"])


def shape(X, w_t, inv, fs_r=np.inf, lam=LAM, cap=True):
    n = X.shape[0]; L = int(round(H.WIN / DT)); step = int(round(H.HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(H.NF, DT)
    live = np.abs(X).sum(axis=0) > 0; Wo = np.abs(np.fft.rfft(w_t, H.NF)); Wo = Wo / Wo.max(); band = (f >= H.LO) & (f <= H.HI)
    e = fs_r if cap else np.inf; b = np.clip((e - f) / TAPER, 0, 1); b = 0.5 - 0.5 * np.cos(np.pi * b)       # 1 below e - 5 Hz, 0 above e
    out = np.zeros_like(X); ws = np.zeros(n); W = 1.0 / (1.0 + inv)
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; Fx = np.fft.rfft(X[i0:i1] * hw[:, None], H.NF, axis=0)
        A = np.sqrt(lib.smooth(np.mean(np.abs(Fx[:, live]) ** 2, axis=1), H.NSM)); A = A / A[band].max()
        T = Wo * np.sqrt(np.interp(f, F10, (W[i0:i1] * hw[:, None]).sum(0) / hw.sum())); T = T / T[band].max()
        G0 = T * A / (A ** 2 + SEL["B3"]["eps"] * np.max(A ** 2)); Gb = G0 ** lam
        G = np.where(f > H.HI, G0, np.where(f < H.LO, G0, b * Gb + (1 - b)))
        out[i0:i1] += np.fft.irfft(Fx * G[:, None], H.NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]
