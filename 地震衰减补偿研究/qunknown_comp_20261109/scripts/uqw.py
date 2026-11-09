"""UQ-W: UQ-IQ with a data-driven noise term (PLAN_STAGE3_ZH.md; DESIGN_ZH.md 2.2 original formula).

UQ-IQ (uqiq.py) minimises E|HG - 1|^2 + s E|H - P|^2 with a constant s (70 dB). The Wiener form of the design,
H = S E[G*] / (S E|G|^2 + N), needs N/S, the noise power over the source-level signal power. The coherence SNR of
the recorded data is snr = S |G|^2 / N, so N/S = E|G|^2 / snr. With rho = s + E|G|^2 / snr(t, f):
    H = (E[G*] + gamma rho E[P]) / (E|G|^2 + rho)
gamma = 1: the UQ-IQ form (phase-only correction where noise dominates); gamma = 0: Wiener (suppressed there).
1/snr = 0 and gamma = 1 give uqiq.uqiq exactly. snr(t, f): lib.coherence_snr of the input data in 0.4 s windows
every 0.1 s (the windows of tvw / VDB), interpolated to each output row and to the operator's frequency grid.

usage: uqw.py smoke     equivalence and code-path checks only (no scores)
"""
import sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqiq import np, lib, semi, lib12, TW, G_IQ, GH, operator_rows, load  # noqa: E402

DT = float(lib12.DT); WIN, HOP, SNR_MIN = 0.4, 0.1, 1e-3


def snr_rows(D, tw, f_op):
    """1/snr at every row time and operator frequency (n, nf)."""
    starts = np.arange(tw[0], tw[-1] - WIN + 1e-9, HOP); cen, inv = [], []
    for a in starts:
        f, snr, _ = lib.coherence_snr(D, tw, a, a + WIN); cen.append(a + WIN / 2)
        inv.append(np.interp(f_op, f, 1 / np.maximum(snr, SNR_MIN)))
    inv = np.array(inv); cen = np.array(cen)
    return np.stack([np.array([np.interp(t, cen, inv[:, k]) for k in range(len(f_op))]) for t in tw])


def expectations(x_hat, dx_hat, sig, f, lf, law):
    """A = E[G*] / E|G|^2, C = 1 / E|G|^2, P = E[P] (|P| <= 1); log domain, as uqiq.operator_rows. (n, nf) each."""
    k = 2 * f * lf; x_hat = np.asarray(x_hat, float)[:, None]; dx_hat = np.asarray(dx_hat, float)[:, None]
    if law == "gauss":
        v = (sig * dx_hat) ** 2; bc = -np.pi * f - 1j * k
        lG2 = -2 * np.pi * f * x_hat + 2 * np.pi ** 2 * f ** 2 * v
        return np.exp(bc * x_hat + 0.5 * bc ** 2 * v - lG2), np.exp(-lG2), np.exp(-1j * k * x_hat - 0.5 * k ** 2 * v)
    e, w = GH; lw = np.log(w / w.sum()); xs = (x_hat + dx_hat * (np.exp(sig * e) - 1)[None, :])[:, None, :]
    l2 = -2 * np.pi * f[None, :, None] * xs + lw; mx = l2.max(-1, keepdims=True)
    lG2 = mx[..., 0] + np.log(np.exp(l2 - mx).sum(-1))
    A = np.exp((-np.pi * f - 1j * k)[None, :, None] * xs + lw - lG2[..., None]).sum(-1)
    P = np.exp(-1j * k[None, :, None] * xs + lw).sum(-1)
    return A, np.exp(-lG2), P


def rows_w(x_hat, dx_hat, sig, f, lf, s, inv_snr, gamma, law="lognormal"):
    """H = (E[G*] + gamma rho E[P]) / (E|G|^2 + rho), rho = s + E|G|^2 inv_snr; divided through by E|G|^2."""
    A, C, P = expectations(x_hat, dx_hat, sig, f, lf, law)
    return (A + gamma * (s * P * C + P * inv_snr)) / (1 + s * C + inv_snr)


def uqw(D, cm, Q, c, sig, tw, t_anchor, gamma, law="lognormal", g_db=G_IQ, nfft=1024, inv_snr=None):
    Q = np.asarray(Q, float); n = len(tw); f = np.fft.rfftfreq(nfft, DT); lag = np.arange(n)
    x_hat = lib.absorption(dict(cm, Q=Q * c), tw); s = semi.gain_sigma2(g_db)
    dx_hat = x_hat - lib.absorption(dict(cm, Q=Q * c), [t_anchor])[0]
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    if inv_snr is None:
        inv_snr = snr_rows(D, tw, f)
    H = rows_w(x_hat, dx_hat, sig, f, lf, s, inv_snr, gamma, law); M = np.empty((n, n))
    for i in range(n):
        h = H[i].copy(); h[0] = h[0].real; h[-1] = h[-1].real
        M[i] = np.fft.irfft(h, nfft)[(i - lag) % nfft]
    return M @ D


if __name__ == "__main__" and sys.argv[1:] == ["smoke"]:
    import uqiq
    from fdcheck import load_fd, uq_t
    t0 = time.monotonic(); arr, meta = load("ZA_c0.75_r0"); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); D = arr["D"]
    f = np.fft.rfftfreq(1024, DT); zero = np.zeros((len(TW), len(f)))
    for law in ("lognormal", "gauss"):
        for sig in (0.0, 0.447):
            a = uqiq.uqiq(D, cm, Q, 0.8, sig, law=law); b = uqw(D, cm, Q, 0.8, sig, TW, 1.4, 1.0, law=law, inv_snr=zero)
            print("semi %-9s sig %.3f  uqw(1/snr=0, gamma=1) vs uqiq rel %.1e" % (law, sig, np.abs(a - b).max() / np.abs(a).max()))
    inv = snr_rows(D, TW, f); fb = (f >= 10) & (f <= 60)
    print("semi ZA_c0.75_r0 coherence SNR 10-60 Hz median: 1.2-1.6 s %.2f, 2.3-2.8 s %.2f" %
          (np.median(1 / inv[(TW >= 1.2) & (TW <= 1.6)][:, fb]), np.median(1 / inv[(TW >= 2.3) & (TW <= 2.8)][:, fb])))
    D0, T, cmf, tw, fp = load_fd("H1"); zf = np.zeros((len(tw), len(f)))
    a = uq_t(D0, cmf, cmf["Q"], 1.0, 0.447, tw, 0.0); b = uqw(D0, cmf, cmf["Q"], 1.0, 0.447, tw, 0.0, 1.0, inv_snr=zf)
    print("FD H1 uqw(1/snr=0, gamma=1) vs uq_t rel %.1e" % (np.abs(a - b).max() / np.abs(a).max()))
    for g in (0.0, 1.0):
        x = uqw(D, cm, Q, 0.8, 0.447, TW, 1.4, g); print("gamma %g finite %s" % (g, bool(np.isfinite(x).all())))
    print("%.0fs" % (time.monotonic() - t0))
