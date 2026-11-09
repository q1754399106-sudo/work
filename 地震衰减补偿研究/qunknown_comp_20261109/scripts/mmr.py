"""Minimax-regret compensation over whole Q-profile scenarios (PLAN_MMR_ZH.md).

Loss of operator row value H at cell c = (t, f) when scenario j is true (the UQ-W noise model, gamma = 0):
    l_cj(H) = |H G_cj - 1|^2 + rho_cj |H|^2,   rho_cj = s + |G_cj|^2 / snr(t, f),
    G_cj = exp[(-pi f + i k) x_j(t)],  k = 2 f ln(f / f_ref)  (G* = exp[(-pi f - i k) x], the uqw.py convention),
    l*_cj = min_H l_cj = rho_cj / (|G_cj|^2 + rho_cj).
Scenarios are whole attenuation profiles after the anchor (where the wavelet is known):
    scalar family  x_j(t) = x_hat(t) + dx_hat(t) (u_j - 1),           u_j on 17 log-spaced points in [0.5, 2];
    two-segment    the multiplier is u_a above t_mid and u_b below it (9 x 9 points), t_mid where dx_hat reaches
                   half its final value; this keeps x_j continuous and non-decreasing.
Design rules (all choose a prior pi over the scenarios; for a fixed pi the minimiser of sum_j pi_j l_cj is
    H_c = sum_j pi_j G*_cj / sum_j pi_j (|G_cj|^2 + rho_cj)                      (Bayes form)
and by Sion's theorem each minimax problem below is solved by this form at its least-favourable pi):
    UQW  pi = lognormal Gauss-Hermite weights (sigma 0.447); identical to uqw.uqw gamma 0 (checked in smoke)
    MMJ  one pi for the section, max_pi sum_c min_H sum_j pi_j l_cj                 (joint worst case)
    RGC  one pi per cell,          max_pi_c min_H sum_j pi_j (l_cj - l*_cj)        (cell-wise regret)
    RGJ  one pi for the section, max_pi sum_c min_H sum_j pi_j (l_cj - l*_cj)      (joint regret, scalar family)
    RGJ2 as RGJ on the two-segment family (pi from every sixth row)
The dual objective phi(pi) is concave; its gradient in pi_j is the (regret) loss of scenario j at H(pi) (Danskin),
so it is maximised by exponentiated-gradient ascent with backtracking. Duality gap max_j R_j - sum_j pi_j R_j is
reported (R_j = mean over cells of l_cj - l*_cj at H(pi)).

usage: mmr.py smoke    equivalence, convergence and timing checks only (no scores)
"""
import sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqw import np, lib, semi, DT, G_IQ, GH, snr_rows, uqw  # noqa: E402

U1 = np.exp(np.linspace(np.log(0.5), np.log(2.0), 17)); U2 = np.exp(np.linspace(np.log(0.5), np.log(2.0), 9))
SIG = 0.447; NFFT = 1024


def setup(cm, Q, c, tw, t_anchor, nfft=NFFT):
    Q = np.asarray(Q, float); f = np.fft.rfftfreq(nfft, DT)
    x_hat = lib.absorption(dict(cm, Q=Q * c), tw); dx_hat = x_hat - lib.absorption(dict(cm, Q=Q * c), [t_anchor])[0]
    lf = np.log(np.maximum(f, np.finfo(float).tiny) / float(cm["fref"])); return f, lf, x_hat, dx_hat


def profiles(x_hat, dx_hat, family):
    """(n, J) attenuation profiles and the (J, 2) multipliers."""
    if family == "scalar":
        return x_hat[:, None] + dx_hat[:, None] * (U1 - 1)[None, :], np.c_[U1, U1]
    half = 0.5 * dx_hat[-1]; i_mid = int(np.argmax(dx_hat >= half)); ua, ub = [a.ravel() for a in np.meshgrid(U2, U2, indexing="ij")]
    lo = np.minimum(dx_hat, dx_hat[i_mid])[:, None]; hi = (dx_hat - np.minimum(dx_hat, dx_hat[i_mid]))[:, None]
    return x_hat[:, None] + lo * (ua - 1)[None, :] + hi * (ub - 1)[None, :], np.c_[ua, ub]


def cells(xs, f, lf, s, inv):
    """G* (n, nf, J), |G|^2 + rho, and l* for every cell and scenario."""
    k = 2 * f * lf; Gc = np.exp((-np.pi * f - 1j * k)[None, :, None] * xs[:, None, :])
    a = np.exp(-2 * np.pi * f[None, :, None] * xs[:, None, :]); rho = s + a * inv[..., None]
    return Gc, a + rho, rho / (a + rho)


def h_of(pi, Gc, den):
    """Bayes-form operator for pi (J,) shared or (n, nf, J) per cell."""
    return (Gc * pi).sum(-1) / (den * pi).sum(-1)


def losses(H, Gc, den, lstar):
    """l_cj(H) - l*_cj (n, nf, J); l = |H|^2 (|G|^2 + rho) - 2 Re(H G) + 1 with G = conj(G*)."""
    return (np.abs(H)[..., None] ** 2 * den - 2 * (H[..., None] * np.conj(Gc)).real + 1) - lstar


def solve_joint(Gc, den, lstar, its=400, tol=1e-6):
    """max over one pi of mean_c min_H sum_j pi_j (l_cj - lstar_cj); lstar = 0 gives the worst case."""
    J = Gc.shape[-1]; pi = np.full(J, 1.0 / J); eta = 10.0

    def ev(p):
        R = losses(h_of(p, Gc, den), Gc, den, lstar).mean((0, 1)); return R, float(p @ R)
    R, phi = ev(pi)
    for it in range(its):
        q = pi * np.exp(eta * (R - R.max())); q /= q.sum(); Rq, pq = ev(q)
        if pq >= phi - 1e-15:
            pi, R, phi = q, Rq, pq; eta *= 1.5
        else:
            eta *= 0.3
        if R.max() - phi < tol * max(abs(phi), 1e-12) or eta < 1e-8:
            break
    return pi, dict(gap=float(R.max() - phi), value=float(phi), its=it + 1, regret=R.tolist())


def solve_cells(Gc, den, lstar, its=300, block=64):
    """max over pi_c per cell, in frequency blocks (memory)."""
    out = [_solve_cells(Gc[:, a:a + block], den[:, a:a + block], lstar[:, a:a + block], its) for a in range(0, Gc.shape[1], block)]
    pi = np.concatenate([o[0] for o in out], 1); gap = np.concatenate([o[1] for o in out], 1); phi = np.concatenate([o[2] for o in out], 1)
    return pi, dict(gap_median=float(np.median(gap)), gap_p99=float(np.quantile(gap, 0.99)),
                    gap_rel_p99=float(np.quantile(gap / np.maximum(phi, 1e-12), 0.99)))


def _solve_cells(Gc, den, lstar, its):
    """vectorised exponentiated gradient with per-cell backtracking."""
    n, nf, J = Gc.shape; pi = np.full((n, nf, J), 1.0 / J); eta = np.full((n, nf, 1), 10.0)
    R = losses(h_of(pi, Gc, den), Gc, den, lstar); phi = (pi * R).sum(-1, keepdims=True)
    for _ in range(its):
        q = pi * np.exp(eta * (R - R.max(-1, keepdims=True))); q /= q.sum(-1, keepdims=True)
        Rq = losses(h_of(q, Gc, den), Gc, den, lstar); pq = (q * Rq).sum(-1, keepdims=True); ok = pq >= phi - 1e-15
        pi = np.where(ok, q, pi); R = np.where(ok, Rq, R); phi = np.where(ok, pq, phi); eta = np.where(ok, eta * 1.5, eta * 0.3)
    return pi, R.max(-1) - phi[..., 0], phi[..., 0]


def apply_rows(H, D, nfft=NFFT):
    n = H.shape[0]; lag = np.arange(n); M = np.empty((n, n))
    for i in range(n):
        h = H[i].copy(); h[0] = h[0].real; h[-1] = h[-1].real
        M[i] = np.fft.irfft(h, nfft)[(i - lag) % nfft]
    return M @ D


def design(D, cm, Q, c, tw, t_anchor, rule, inv=None, g_db=G_IQ):
    """Operator output for one rule; returns (output, info)."""
    f, lf, x_hat, dx_hat = setup(cm, Q, c, tw, t_anchor); s = semi.gain_sigma2(g_db)
    if inv is None:
        inv = snr_rows(D, tw, f)
    if rule == "UQW":
        e, w = GH; xs = x_hat[:, None] + dx_hat[:, None] * (np.exp(SIG * e) - 1)[None, :]
        Gc, den, _ = cells(xs, f, lf, s, inv); return apply_rows(h_of(w / w.sum(), Gc, den), D), {}
    if rule == "RGJ2":
        xs, uu = profiles(x_hat, dx_hat, "two")
        sub = slice(None, None, 6); Gc, den, ls = cells(xs[sub], f, lf, s, inv[sub]); pi, info = solve_joint(Gc, den, ls)
        keep = pi > 0; Gc, den, _ = cells(xs[:, keep], f, lf, s, inv)      # zero-weight profiles drop out exactly (memory)
        info["pi"] = pi.tolist(); info["u"] = uu.tolist(); return apply_rows(h_of(pi[keep], Gc, den), D), info
    xs, uu = profiles(x_hat, dx_hat, "scalar")
    if rule in ("MMJ", "RGJ"):                      # the shared pi from every second row; the operator on all rows
        Gc, den, ls = cells(xs[::2], f, lf, s, inv[::2])
        pi, info = solve_joint(Gc, den, np.zeros_like(ls) if rule == "MMJ" else ls)
    Gc, den, ls = cells(xs, f, lf, s, inv)
    if rule == "RGC":
        pi, info = solve_cells(Gc, den, ls); info["pi_mean"] = pi.mean((0, 1)).tolist()
        return apply_rows(h_of(pi, Gc, den), D), info
    elif rule not in ("MMJ", "RGJ"):
        raise ValueError(rule)
    info["pi"] = pi.tolist(); info["u"] = uu[:, 0].tolist(); return apply_rows(h_of(pi, Gc, den), D), info


if __name__ == "__main__" and sys.argv[1:] == ["smoke"]:
    from uqiq import TW, load
    import fdcheck as fd
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); arr, meta = load("ZA_c0.75_r0"); D = arr["D"]
    inv = snr_rows(D, TW, np.fft.rfftfreq(NFFT, DT))
    for c in (1.0 * meta["c_eff_cum"], 2.0 * meta["c_eff_cum"]):
        a = uqw(D, cm, Q, c, SIG, TW, 1.4, 0.0, inv_snr=inv); b, _ = design(D, cm, Q, c, TW, 1.4, "UQW", inv)
        print("semi c %.3f  Bayes-form UQW vs uqw.uqw rel %.1e" % (c, np.abs(a - b).max() / np.abs(a).max()))
    for rule in ("MMJ", "RGJ", "RGC", "RGJ2"):
        t0 = time.monotonic(); X, info = design(D, cm, Q, meta["c_eff_cum"], TW, 1.4, rule, inv)
        short = {k: (np.round(v, 3).tolist() if isinstance(v, list) and k in ("pi", "pi_mean") else v) for k, v in info.items() if k not in ("regret", "u")}
        print(rule, "%.0fs" % (time.monotonic() - t0), "finite", bool(np.isfinite(X).all()), short)
    D0, T, cmf, tw, fp = fd.load_fd("H1"); Df = D0 + fd.noise(D0, 14.0, 1); invf = snr_rows(Df, tw, np.fft.rfftfreq(NFFT, DT))
    a = uqw(Df, cmf, cmf["Q"], 1.0, SIG, tw, 0.0, 0.0, inv_snr=invf); b, _ = design(Df, cmf, cmf["Q"], 1.0, tw, 0.0, "UQW", invf)
    print("FD H1 Bayes-form UQW vs uqw.uqw rel %.1e" % (np.abs(a - b).max() / np.abs(a).max()))
    for rule in ("RGJ", "MMJ"):
        X, info = design(Df, cmf, cmf["Q"], 1.0, tw, 0.0, rule, invf); print("FD", rule, np.round(info["pi"], 3).tolist(), "gap %.1e" % info["gap"])
