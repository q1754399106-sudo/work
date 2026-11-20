"""Shared code of anchor_active_hessian_20261120 (PLAN_ZH.md, DERIVATION_ZH.md).

The QScan lasso with column-normalised dictionary B = B~/scale and weights pen_k = lam sqrt(2 ln n) sigma
median(scale)/scale_k is the same problem as an unweighted L1 on the unnormalised dictionary B~ = G Win with the
single weight mu sigma, mu = lam sqrt(2 ln n) median(scale) (y = x/scale). On a fixed active set A with signs s the
per-trace value is F = 1/2 |d|^2 - 1/2 c' (B~_A' B~_A)^-1 c, c = B~_A' d - mu sigma s, and the section objective is
J = sum_j F_j / sigma_j^2 (live traces), identical to the scan's J.

theta_m = extra cumulative t* of the propagation column m of G (G[:, m] = irfft(exp(b (a_m v + theta_m)) shift_m)), so
d/dv = sum_m a_m d/dtheta_m at fixed mu. Analytic Hessian of J in theta on the fixed active set (mu fixed):
    H = sum_j w_j [ Z (G1'G1) Z - diag(z o G2'r) - Hc' (B~_A'B~_A)^-1 Hc ],
    Hc = (B~_A' G1) Z - Win[:, A]' diag(G1' r),   z = Win[:, A] y_A,  r = d - B~_A y_A,
G1, G2 the first and second theta-derivatives of the columns of G. The first two terms are the "direct" part (fixed
coefficients, active energy only); the last is the coupling (re-fit of the active coefficients).
"""
import os, sys, json
from pathlib import Path

for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

HERE = Path(__file__).resolve().parents[1]
RESEARCH = HERE.parent
_CACHE = Path("field_baselines_volve_20261029") / "cache" / "inputs"
DATA_ROOT = RESEARCH if (RESEARCH / _CACHE).exists() else Path(os.environ.get("QRES_ROOT", r"D:\ChatGPT\研究\地震衰减补偿研究"))
ANCH_CODE = DATA_ROOT / "qscan_anchor_theory_20261030" / "scripts"     # env15 resolves the input caches from here
ANCH = RESEARCH / "qscan_anchor_theory_20261030"                         # archived scan results (read only)
sys.path.insert(0, str(ANCH_CODE))
import qprofile as qp  # noqa: E402
from env15 import np, lib, semi, sc11, bm, gen, E12OUT  # noqa: E402

LAM = 0.5
PROFILES = {"S": [30.0, 200.0, 200.0, 200.0], "D": [250.0, 250.0, 40.0, 40.0],
            "D2": [150.0, 150.0, 60.0, 60.0], "M": [200.0, 45.0, 200.0, 200.0],
            "R1": [100.0, 100.0, 100.0, 100.0], "R2": [60.0, 60.0, 150.0, 150.0], "R3": [150.0, 40.0, 150.0, 150.0],
            "R4": [50.0, 90.0, 250.0, 250.0], "X1": [120.0, 120.0, 70.0, 70.0], "X2": [300.0, 300.0, 60.0, 60.0],
            "X3": [150.0, 150.0, 150.0, 50.0]}
SEEDS = {"S": (9971, 9972, 9973, 9974), "D": (9971, 9972, 9973, 9974), "D2": (9981, 9982, 9983, 9984),
         "M": (9981, 9982, 9983, 9984)}
SCANDIR = {"S": "results7", "D": "results7", "D2": "results10", "M": "results10"}
qp.PROFILES.clear(); qp.PROFILES.update(PROFILES)


def datasets():
    out = [("base", f"T_c{c}_s{s}", c, s) for c in (0.65, 0.85, 1.05, 1.25, 1.5) for s in (9701, 9702)]
    for p in PROFILES:
        for s in SEEDS.get(p, (9991, 9992)):
            for c in qp.CS:
                out.append((p, f"c{c}_s{s}", c, s))
    return out


def load(prof, did, c, s, base):
    """Section D and model cm exactly as the archived scans (anchor.py / qprofile.py / bound.py)."""
    if prof == "base":
        return np.load(E12OUT / f"{did}.npz")["Input"].astype(float), semi.curtain_cm()
    ctx = qp.context(base, prof); g = gen.signal(ctx, c, s)
    return g["x"] + qp.MULT * g["level"] * ctx["noise"][0], ctx["cm"]


def v_selected(prof, did):
    """Selection at anchor c_a = 1 from the archived scan: parabola vertex in v through the grid minimum."""
    VG = 1.0 / np.asarray(lib.QGRID, float)
    if prof == "base":
        J = np.asarray(json.loads((ANCH / "results" / f"{did}__a1.0.json").read_text(encoding="utf-8"))["diag"]["J"][0], float)
        jm = int(np.argmin(J)); v = float(VG[jm])
        if 0 < jm < len(J) - 1:
            idx = [jm - 1, jm, jm + 1]; p2 = np.polyfit(VG[idx], J[idx], 2)
            if p2[0] > 0:
                v = float(np.clip(-p2[1] / (2 * p2[0]), VG[jm + 1], VG[jm - 1]))
        return v
    d = SCANDIR.get(prof, "results26")
    return float(json.loads((ANCH / d / f"{prof}_{did}__a1.0.json").read_text(encoding="utf-8"))["v_sel"])


def freq(cm):
    N = int(cm["nfft"]); f = np.arange(N // 2 + 1) / (N * float(cm["dt"]))
    b = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); b[0] = 0
    return N, f, b


def gmats(cm, aeff, order=0, aamp=None):
    """Columns G[:, m] = irfft(exp(b aeff_m) shift_m) and (order >= 1, 2) their first/second aeff-derivatives.
    aamp (order 0 only): separate t* for the amplitude part -pi f of b."""
    t = bm.vec(cm["t"]); n = len(t); N, f, b = freq(cm)
    if aamp is None:
        E0 = np.exp(np.outer(b, aeff))
    else:
        E0 = np.exp(np.outer(b.real, aamp) + 1j * np.outer(b.imag, aeff))
    E = E0 * np.exp(-2j * np.pi * np.outer(np.arange(N // 2 + 1), np.arange(n)) / N)
    out = []
    for p in range(order + 1):
        M = E * (b[:, None] ** p); M[-1] = M[-1].real
        out.append(np.fft.irfft(M, n=N, axis=0)[:n])
    return out


def wavelet(D, cm, va):
    """Anchor wavelet: lib.wavelets with the assumed Q = model Q / va (cumulative t* at 1.4 s = a0 va)."""
    Q = np.asarray(cm["Q"], float)
    return lib.wavelets(D, semi.TW, dict(cm, Q=Q / va))["w_in"]


def wavelet_raw(D, cm, va, nfft=512):
    """lib.wavelets w_in without the division by max|w| (the global factor then follows the anchor physically).
    With mu recomputed from the dictionary the scan objective is invariant to this factor."""
    Q = np.asarray(cm["Q"], float); c = dict(cm, Q=Q / va)
    f, snr_s, P_s = lib.coherence_snr(D, semi.TW, 1.20, 1.60, nfft)
    _, snr_d, _ = lib.coherence_snr(D, semi.TW, 2.00, 3.00, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * lib.absorption(c, [1.40])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))]); fs = lib.snr_band(f, snr_s, fpk)
    g = np.ones_like(f); r = (f > fs - 5) & (f < fs)
    g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fs - 5)) / 5)); g[f >= fs] = 0; g[f < 3] = 0
    return np.roll(np.fft.irfft(amp * g, nfft), 100)[:201]


def unnorm_dict(cm, v, w_in, theta=None, theta_amp=None):
    """B~ = G Win; theta adds t* to amplitude and phase of the propagation columns, theta_amp to the amplitude only."""
    a = lib.absorption(cm, bm.vec(cm["t"])); aeff = a * v + (0 if theta is None else theta)
    G = gmats(cm, aeff, aamp=None if theta_amp is None else aeff + theta_amp)[0]; Win = lib.toeplitz_wavelet(w_in, len(a))
    return G @ Win, Win


def mu_of(Bt, n, q=50.0):
    """L1 weight of the scan: lam sqrt(2 ln n) x median column norm (q = 50; other percentiles for the intervention)."""
    s = np.linalg.norm(Bt, axis=0)
    return LAM * np.sqrt(2 * np.log(n)) * float(np.median(s) if q == 50.0 else np.percentile(s, q))


def solve(D, Bt, pp, start=None, tol=None):
    """Pipeline lasso on the normalised dictionary; returns y = X/scale, J (scan objective), converged."""
    n = D.shape[0]; scale = np.linalg.norm(Bt, axis=0); B = Bt / scale
    old = bm.TOL
    if tol is not None:
        bm.TOL = tol
    try:
        X, info = bm.solve(B, np.eye(n), D, pp["sigma"], scale, pp["K"], LAM, 0.0, pp["live"],
                           start=None if start is None else start * scale[:, None])
    finally:
        bm.TOL = old
    sigma, live = pp["sigma"], pp["live"]; s2 = np.where(live, sigma ** 2, np.inf)
    pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]; res = D - B @ X
    J = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + LAM * np.sum(pen1 * np.abs(X), axis=0))[live] / s2[live]))
    return X / scale[:, None], J, dict(converged=bool(info["converged"]), iterations=int(info["iterations"]),
                                       kkt=float(info["kkt_ratio"]))


def polish(D, Bt, y, pp, mu):
    """Exact fixed-active-set solution on the support/signs of y; reports sign and KKT consistency."""
    live = np.where(pp["live"])[0]; sigma = pp["sigma"]; A, S, Y = {}, {}, np.zeros_like(y)
    bad_sign = 0; nact = 0; kkt = []
    for j in live:
        Aj = np.flatnonzero(y[:, j] != 0); s = np.sign(y[Aj, j]); BA = Bt[:, Aj]
        yA = np.linalg.solve(BA.T @ BA, BA.T @ D[:, j] - mu * sigma[j] * s) if len(Aj) else np.zeros(0)
        bad_sign += int(np.sum(np.sign(yA) != s)); nact += len(Aj)
        Y[Aj, j] = yA; A[j], S[j] = Aj, s
        r = D[:, j] - BA @ yA; g = np.abs(Bt.T @ r) / (mu * sigma[j]); g[Aj] = 0; kkt.append(float(g.max()))
    return Y, A, S, dict(n_active=nact, bad_sign=bad_sign, kkt_inactive_max=float(np.max(kkt)),
                         kkt_inactive_p99=float(np.percentile(kkt, 99)))


def jfix(D, Bt, A, S, pp, mu):
    """Fixed-active-set closed-form objective J (same scaling as the scan objective)."""
    J = 0.0
    for j, Aj in A.items():
        d = D[:, j]; sg = pp["sigma"][j]
        if len(Aj):
            BA = Bt[:, Aj]; c = BA.T @ d - mu * sg * S[j]
            F = 0.5 * d @ d - 0.5 * c @ np.linalg.solve(BA.T @ BA, c)
        else:
            F = 0.5 * d @ d
        J += F / sg ** 2
    return J


def theta_hessian(D, cm, v, w_in, Y, A, pp):
    """Analytic Hessian of J in theta (mu fixed) on the fixed active set: total H, direct part Hd; also the wavefield
    energy profile z^2 (weighted like J)."""
    a = lib.absorption(cm, bm.vec(cm["t"])); n = len(a)
    G, G1, G2 = gmats(cm, a * v, order=2); Win = lib.toeplitz_wavelet(w_in, n); Bt = G @ Win
    M11 = G1.T @ G1; BG1 = Bt.T @ G1
    H = np.zeros((n, n)); Hd = np.zeros((n, n)); z2 = np.zeros(n)
    for j, Aj in A.items():
        if not len(Aj):
            continue
        w = 1.0 / pp["sigma"][j] ** 2; yA = Y[Aj, j]; z = Win[:, Aj] @ yA; r = D[:, j] - Bt[:, Aj] @ yA
        g1r = G1.T @ r; g2r = G2.T @ r
        direct = M11 * np.outer(z, z) - np.diag(z * g2r)
        Hc = BG1[Aj, :] * z[None, :] - Win[:, Aj].T * g1r[None, :]
        BA = Bt[:, Aj]; coup = Hc.T @ np.linalg.solve(BA.T @ BA, Hc)
        Hd += w * direct; H += w * (direct - coup); z2 += w * z ** 2
    return H, Hd, z2


def alpha_full(H, a, a0):
    return float(a0 * (a @ H @ np.ones_like(a)) / (a @ H @ a))


def alpha_diag(phi, a, a0):
    phi = np.asarray(phi, float); return float(a0 * np.sum(phi * a) / np.sum(phi * a ** 2))


def smooth(x, w=25):
    x = np.asarray(x, float); k = np.ones(w)
    return np.convolve(x, k, "same") / np.convolve(np.ones_like(x), k, "same")
