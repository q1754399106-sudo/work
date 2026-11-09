"""UQ-IQ: inverse Q under an uncertain Q factor (DESIGN_ZH.md section 2, PLAN_STAGE1_ZH.md).

The 70 dB stabilised inverse Q of experiment 12 (lib12.iq with c_amp = c_ph = c),
    H = exp(-i k x) (beta + s) / (beta^2 + s),  beta = exp(-pi f x),  k = 2 f ln(f / f_ref),  x = t*(t) at c,
is the minimiser of |H G - 1|^2 + s |H - P|^2 with G = exp[(-pi f + i k) x] (the forward model) and P = exp(-i k x)
(phase-only correction). With x uncertain, minimising the expectation gives
    H = (E[G*] + s E[P]) / (E|G|^2 + s),
which is exactly lib12.iq when the uncertainty is zero.

Only attenuation relative to the wavelet window is uncertain. The wavelet is de-attenuated with the model t* at
1.4 s (lib.wavelets), and the 1.4 s spectrum of the data does not depend on the true Q, so after compensation with
u_hat = 1/c the output spectrum differs from the wavelet by exp{-pi f [da(t) (u - u_hat) + a(1.4) (1 - u_hat)]},
da = a(t) - a(1.4). The second term is fixed by u_hat; only the first is uncertain. Hence
    x = x_hat + dx_hat (u / u_hat - 1),  dx_hat = x_hat(t) - x_hat(1.4)   (x_hat = t* at c, as lib12.iq)
with relative uncertainty sig of u (log scale).

law "gauss": u / u_hat - 1 ~ N(0, sig^2), v = (sig dx_hat)^2, closed form (Gaussian moment generating function):
    E[G*] = exp(b* x_hat + b*^2 v / 2),  E|G|^2 = exp(-2 pi f x_hat + 2 pi^2 f^2 v),  E[P] = exp(-i k x_hat - k^2 v / 2)
law "lognormal": u / u_hat = exp(sig e), e ~ N(0, 1), 41-point Gauss-Hermite quadrature (log-sum-exp).
law "gauss_abs": as gauss but v = (sig x_hat)^2 (uncertainty on the absolute t*; mis-specified here, descriptive).
Evaluated as (E[G*] / E|G|^2 + s E[P] / E|G|^2) / (1 + s / E|G|^2) to avoid overflow.

usage: uqiq.py smoke      sig = 0 must reproduce lib12.iq (max |d| / max |IQ70| <= 1e-10) on two datasets x 3 m
"""
import sys, json, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "mu_tune_20261031" / "scripts"))
import lib17  # noqa: E402
from lib17 import np, lib, semi, lib12, TW, G_IQ, load  # noqa: E402

DT = float(lib12.DT)
GH = np.polynomial.hermite_e.hermegauss(41)
T_ANCHOR = 1.40                                     # lib.wavelets de-attenuates the wavelet with the model t* here


def operator_rows(x_hat, dx_hat, sig, f, lf, s, law="gauss"):
    """H for each time row (x_hat: t* per row at c; dx_hat: x_hat minus its value at the anchor). (n, nf) complex."""
    k = 2 * f * lf; x_hat = np.asarray(x_hat, float)[:, None]; dx_hat = np.asarray(dx_hat, float)[:, None]
    if law in ("gauss", "gauss_abs"):
        v = (sig * (dx_hat if law == "gauss" else x_hat)) ** 2; bc = -np.pi * f - 1j * k
        lG2 = -2 * np.pi * f * x_hat + 2 * np.pi ** 2 * f ** 2 * v
        lGc = bc * x_hat + 0.5 * bc ** 2 * v
        lP = -1j * k * x_hat - 0.5 * k ** 2 * v
        A = np.exp(lGc - lG2); B = np.exp(lP - lG2); C = np.exp(-lG2)
    elif law == "lognormal":
        # log domain: in the shallow rows (da < 0) a large sig drives the tail nodes to large negative x
        e, w = GH; lw = np.log(w / w.sum()); xs = (x_hat + dx_hat * (np.exp(sig * e) - 1)[None, :])[:, None, :]  # (n,1,q)
        l2 = -2 * np.pi * f[None, :, None] * xs + lw; mx = l2.max(-1, keepdims=True)
        L2 = mx[..., 0] + np.log(np.exp(l2 - mx).sum(-1))                                       # log E|G|^2
        A = np.exp((-np.pi * f - 1j * k)[None, :, None] * xs + lw - L2[..., None]).sum(-1)       # E[G*] / E|G|^2
        B = np.exp(-1j * k[None, :, None] * xs + lw - L2[..., None]).sum(-1)                    # E[P] / E|G|^2
        C = np.exp(-L2)
    else:
        raise ValueError(law)
    return (A + s * B) / (1 + s * C)


def uqiq(D, cm, Q, c, sig, g_db=G_IQ, law="gauss", nfft=1024):
    """Same structure as lib12.iq (row i: filter for time TW[i], applied as a time-varying convolution)."""
    Q = np.asarray(Q, float); n = len(TW); f = np.fft.rfftfreq(nfft, DT); lag = np.arange(n)
    x_hat = lib.absorption(dict(cm, Q=Q * c), TW); s = semi.gain_sigma2(g_db)
    dx_hat = x_hat - lib.absorption(dict(cm, Q=Q * c), [T_ANCHOR])[0]
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    H = operator_rows(x_hat, dx_hat, sig, f, lf, s, law); M = np.empty((n, n))
    for i in range(n):
        h = H[i].copy(); h[0] = h[0].real; h[-1] = h[-1].real
        M[i] = np.fft.irfft(h, nfft)[(i - lag) % nfft]
    return M @ D


def gain_db(cm, Q, c, sig, t, law="gauss", nfft=1024):
    """|H| in dB at time t (descriptive: the uncertainty high-cut)."""
    f = np.fft.rfftfreq(nfft, DT); lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"]))
    x = lib.absorption(dict(cm, Q=np.asarray(Q, float) * c), [t, T_ANCHOR])
    return f, 20 * np.log10(np.abs(operator_rows(x[:1], x[:1] - x[1], sig, f, lf, semi.gain_sigma2(G_IQ), law)[0]) + 1e-300)


if __name__ == "__main__" and sys.argv[1:] == ["smoke"]:
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); worst = 0.0; t0 = time.monotonic()
    for did in ("ZA_c0.75_r0", "ZF_c1.05_r1"):
        arr, meta = load(did)
        for m in (0.5, 1.0, 2.0):
            c = m * meta["c_eff_cum"]; ref = lib12.iq(arr["D"], cm, Q, c, c, G_IQ)
            for law in ("gauss", "lognormal", "gauss_abs"):
                r = float(np.abs(uqiq(arr["D"], cm, Q, c, 0.0, law=law) - ref).max() / np.abs(ref).max())
                worst = max(worst, r); print(did, m, law, "rel max|d| %.2e" % r)
    # descriptive operator property (no scores): first frequency where |H| at 2.55 s is 6 / 20 dB below IQ70
    first = lambda f, ok: "%.1f Hz" % f[np.argmax(ok)] if ok.any() else "none"
    for c in (0.75, 1.05):
        f, g0 = gain_db(cm, Q, c, 0.0, 2.55)
        for law in ("gauss", "lognormal", "gauss_abs"):
            for sig in (0.1225, 0.245, 0.447):
                d = gain_db(cm, Q, c, sig, 2.55, law)[1] - g0
                print("c %.2f %-9s sig %.3f: -6 dB at %s, -20 dB at %s" % (c, law, sig, first(f, (d < -6) & (f > 5)),
                                                                          first(f, (d < -20) & (f > 5))))
    print("worst %.2e  %.0fs" % (worst, time.monotonic() - t0)); sys.exit(0 if worst <= 1e-10 else 1)
