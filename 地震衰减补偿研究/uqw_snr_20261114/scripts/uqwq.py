"""UQW-Q: UQ-W with a posterior (data-updated) attenuation spectrum, built into the operator (development candidate).

Measurement, per Hann window (0.5 s, hop 0.1 s) centred at t, from the input data only:
  r(f) = 0.5 ln(C_pred / C): C the coherent power P / (1 + 1/snr); C_pred the anchor signal power attenuated with the input Q
  (uqwp.py). Reflectivity colour makes r rough in f, a wrong Q (or a smooth non-constant-Q change) makes it smooth, so only a
  weighted quadratic in f is kept: r_s(f) = a1 f + a2 f^2 (intercept, i.e. the level, discarded), weights w w0 inside LO-HI Hz
  (w = snr / (1 + snr) in the window, w0 at the anchor); held constant outside LO-HI.
  v_d(f): variance of r_s(f) from the residual scatter of the fit (weighted least squares, KCELL bins per independent cell).
Update, no free parameter:
  lam(t, f) = v_p / (v_p + v_d),  v_p = (pi f |dx_hat(t)| sig)^2  (prior variance of the log amplitude; zero at the anchor)
  posterior mean: log amplitude of G shifted by -lam r_s;   posterior spread: sig_eff(t, f) = sig sqrt(1 - lam).
Operator:  H = H_uqw(sig_eff(t, f)) exp(lam r_s)   (uqw.rows_w with a per-cell sig; lam = 0 gives uqw.uqw exactly).
  mode "mean": only the mean shift (sig kept);  mode "post": mean shift and sig_eff.
"""
from env import np, lib, DT, H, CW, F10, semi
from uqw import GH, G_IQ

LO, HI, NF, NSM, KCELL, RMAX = H.LO, H.HI, H.NF, H.NSM, 11.0, 3.0


def anchor_semi(D, tw, a=1.2, b=1.6):
    f, snr, P = lib.coherence_snr(D, tw, a, b, NF); w0 = snr / (1 + snr)
    return np.sqrt(lib.smooth(P * w0, NSM)), w0


def fields(D, tw, dx_hat, inv, A0, w0, sig):
    """(lam, shift) on (rows, F10): shift = lam r_s, the log-amplitude gain to apply."""
    n = D.shape[0]; L = int(round(H.WIN / DT)); step = int(round(H.HOP / DT)); h = np.hanning(L); f = np.fft.rfftfreq(NF, DT)
    live = np.abs(D).sum(axis=0) > 0; band = (f >= LO) & (f <= HI); W = 1.0 / (1.0 + inv); cen, LAM, RS = [], [], []
    fb = f[band]; Xd = np.stack([np.ones_like(fb), fb, fb ** 2], axis=1); fc = np.clip(f, LO, HI); Xa = np.stack([fc, fc ** 2], axis=1)
    for s0 in range(-L // 2, n - L // 2, step):
        i0, i1 = max(s0, 0), min(s0 + L, n); hw = h[i0 - s0:i1 - s0]
        P = lib.smooth(np.mean(np.abs(np.fft.rfft(D[i0:i1][:, live] * hw[:, None], NF, axis=0)) ** 2, axis=1), NSM)
        w = np.interp(f, F10, (W[i0:i1] * hw[:, None]).sum(0) / hw.sum()); dx = float((dx_hat[i0:i1] * hw).sum() / hw.sum())
        C = np.maximum(P * w, 1e-30 * P.max()); r = np.clip(0.5 * np.log(np.maximum(A0 ** 2, 1e-30 * (A0 ** 2).max()) * np.exp(-2 * np.pi * f * dx) / C), -2 * RMAX, 2 * RMAX)
        wg = np.maximum((w * w0)[band], 1e-6); G = Xd.T @ (wg[:, None] * Xd); Gi = np.linalg.pinv(G); a = Gi @ (Xd.T @ (wg * r[band]))
        s2 = float((wg * (r[band] - Xd @ a) ** 2).sum() / wg.sum()) * KCELL * (Xd.T @ (wg[:, None] ** 2 * Xd)).trace() / max(G.trace() * wg.sum() / len(wg), 1e-30) * 0 + \
            float((wg * (r[band] - Xd @ a) ** 2).sum() / wg.sum())
        cov = Gi * s2 * KCELL * wg.sum() / len(wg)                       # WLS covariance with mean-normalised weights, KCELL bins per cell
        rs = np.clip(Xa @ a[1:], -RMAX, RMAX); vd = np.einsum("ij,jk,ik->i", Xa, cov[1:, 1:], Xa)
        vp = (np.pi * f * abs(dx) * sig) ** 2; lam = vp / (vp + vd + 1e-30); cen.append(tw[0] + DT * (s0 + L / 2)); LAM.append(lam); RS.append(lam * rs)
    cen = np.array(cen); o = np.argsort(cen); cen = cen[o]
    it = lambda A: np.stack([np.interp(F10, f, row) for row in np.stack([[np.interp(t, cen, np.array(A)[o][:, k]) for k in range(len(f))] for t in tw])])
    return it(LAM), it(RS)


def rows(x_hat, dx_hat, sig, f, lf, s, inv):
    """uqw.rows_w (lognormal, gamma 0) with sig an (n, nf) array."""
    e, w = GH; lw = np.log(w / w.sum()); k = 2 * f * lf
    xs = x_hat[:, None, None] + dx_hat[:, None, None] * (np.exp(sig[:, :, None] * e[None, None, :]) - 1)
    l2 = -2 * np.pi * f[None, :, None] * xs + lw; mx = l2.max(-1, keepdims=True); lG2 = mx[..., 0] + np.log(np.exp(l2 - mx).sum(-1))
    A = np.exp((-np.pi * f - 1j * k)[None, :, None] * xs + lw - lG2[..., None]).sum(-1)
    return A / (1 + s * np.exp(-lG2) + inv)


def uqwq(D, cm, Q, c, sig, tw, t_anchor, inv, A0, w0, mode="post", nfft=1024, ret=False):
    Q = np.asarray(Q, float); n = len(tw); f = np.fft.rfftfreq(nfft, DT); lag = np.arange(n)
    x_hat = lib.absorption(dict(cm, Q=Q * c), tw); dx_hat = x_hat - lib.absorption(dict(cm, Q=Q * c), [t_anchor])[0]
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); lam, shift = fields(D, tw, dx_hat, inv, A0, w0, sig)
    se = sig * np.sqrt(1 - lam) if mode == "post" else np.full_like(lam, sig)
    Hm = rows(x_hat, dx_hat, se, f, lf, semi.gain_sigma2(G_IQ), inv) * np.exp(shift); M = np.empty((n, n))
    for i in range(n):
        hh = Hm[i].copy(); hh[0] = hh[0].real; hh[-1] = hh[-1].real; M[i] = np.fft.irfft(hh, nfft)[(i - lag) % nfft]
    return (M @ D, lam, shift) if ret else M @ D
