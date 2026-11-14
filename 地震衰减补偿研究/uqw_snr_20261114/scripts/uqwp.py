"""UQW-P: UQ-W with a data update of the amplitude that needs no target wavelet (development candidate).

In each Hann window (0.5 s, hop 0.1 s, as uqwh.half) centred at t:
  C(t, f)   coherent (signal) power of the input, P / (1 + 1/snr) with the coherence snr that UQ-W already uses;
  C_pred    the anchor signal power A0(f)^2 attenuated with the input Q: A0^2 exp(-2 pi f [x_hat(t) - x_hat(anchor)]);
  r(f)      = 0.5 ln(C_pred / C), minus its weighted mean over LO-HI Hz (the level is not corrected).
r is the departure, in log amplitude, of the measured spectral change since the anchor from the constant-Q prediction; it is
zero when Q is right and the reflectivity is white, and pi f da (u - u_hat) when Q is wrong by a factor. The UQ-W output is
multiplied by exp(lam(t, f) r(f)). lam:
  "const"  LAM inside LO-HI, 0 outside
  "rel"    v_p / (v_p + v_c) w(t, f) w0(f),  v_p = (pi f |dx_hat| sig)^2 the prior variance of the log amplitude,
           v_c(f) the variance of r from reflectivity colour (measured, vcolour.py), w = snr / (1 + snr) in the window and
           w0 the same at the anchor: no free parameter.
A0: semi-synthetic, the coherent amplitude spectrum of the 1.2-1.6 s window of the data; FD, the known source spectrum.
"""
from env import np, lib, TW, DT, H, CW, snr_rows, F10

LAM, LO, HI, NF, NSM, RMAX = 0.5, H.LO, H.HI, H.NF, H.NSM, 3.0


def anchor_semi(D, tw, a=1.2, b=1.6):
    """(A0(f), w0(f)) on rfftfreq(NF): coherent amplitude spectrum and Wiener weight of the anchor window."""
    f, snr, P = lib.coherence_snr(D, tw, a, b, NF); w0 = snr / (1 + snr)
    return np.sqrt(lib.smooth(P * w0, NSM)), w0


def update(X, D, tw, dx_hat, inv, A0, w0, sig, mode="rel", vc=None, lam=LAM, ret=False):
    n = X.shape[0]; L = int(round(H.WIN / DT)); step = int(round(H.HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(NF, DT)
    live = np.abs(D).sum(axis=0) > 0; band = (f >= LO) & (f <= HI); out = np.zeros_like(X); ws = np.zeros(n); W = 1.0 / (1.0 + inv); log = []
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]; Fx = np.fft.rfft(X[i0:i1] * hw[:, None], NF, axis=0)
        P = lib.smooth(np.mean(np.abs(np.fft.rfft(D[i0:i1][:, live] * hw[:, None], NF, axis=0)) ** 2, axis=1), NSM)
        w = np.interp(f, F10, (W[i0:i1] * hw[:, None]).sum(0) / hw.sum()); dx = float((dx_hat[i0:i1] * hw).sum() / hw.sum())
        C = np.maximum(P * w, 1e-30 * P.max()); r = 0.5 * np.log(np.maximum(A0 ** 2, 1e-30 * (A0 ** 2).max()) * np.exp(-2 * np.pi * f * dx) / C)
        wg = (w * w0 * np.sqrt(C))[band]; r = np.clip(r - (r[band] * wg).sum() / wg.sum(), -RMAX, RMAX)
        if mode == "const":
            lm = np.where(band, lam, 0.0)
        else:
            vp = (np.pi * f * abs(dx) * sig) ** 2; lm = (vp / (vp + vc)) * w * w0 if mode == "rel" else w * w0
        out[i0:i1] += np.fft.irfft(Fx * np.exp(lm * r)[:, None], NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw; log.append((tw[0] + DT * (s0 + L / 2), r, lm))
    Y = out / np.maximum(ws, 1e-3 * ws.max())[:, None]
    return (Y, log) if ret else Y


def dxhat(cm, Q, c, tw, t_anchor):
    x = lib.absorption(dict(cm, Q=np.asarray(Q, float) * c), tw); return x - lib.absorption(dict(cm, Q=np.asarray(Q, float) * c), [t_anchor])[0]
