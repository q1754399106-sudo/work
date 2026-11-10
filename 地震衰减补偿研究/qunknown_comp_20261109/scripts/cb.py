"""Restricted-Bayes compensation (PLAN_CB_ZH.md, as amended before any score): expected loss under the lognormal
prior truncated to the 17-point profile grid u in [0.5, 2] (weights exp(-ln^2 u / 2 sigma^2), sigma 0.447),
subject to every profile regret staying below a cap:

    min_H  mean_c sum_j p_j l_cj(H)   s.t.  R_j(H) = mean_c [l_cj(H) - l*_cj] <= B_tau   (all 17 profiles j)
    B_tau = R_mm + tau (R_T - R_mm),

R_mm the minimax-regret value (mmr.solve_joint, the RGJ dual value), R_T = max_j R_j at lambda = 0. The Lagrangian
dual gives the Bayes form with prior p + lambda (lambda_j >= 0):
    H = sum_j (p_j + lambda_j) G*_j / sum_j (p_j + lambda_j)(|G_j|^2 + rho_j)
tau = 1: lambda = 0, the truncated-prior UQ-W (UQWt). The dual g(lambda) = F0(H) + sum_j lambda_j (R_j(H) - B) is
concave with gradient R_j(H(lambda)) - B; projected gradient ascent with backtracking.
Amendment reason: the untruncated Gauss-Hermite prior has nodes up to u = 53; above the anchor (dx_hat < 0) they
are negative-attenuation profiles (x down to -1.25 s) whose |G|^2 (up to e^982) dominates the expected loss despite
weights of 1e-30, so the cap could not be met (smoke: 1.7 % over at tau 0.5, m 1; 3-12 % at tau 0.25).
Rows: every second row for the cap and the solve (as RGJ), all rows (in blocks) for the operator.

usage: cb.py smoke    constraint and timing checks only (no scores)
"""
import sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mmr import np, semi, DT, NFFT, SIG, G_IQ, U1, setup, profiles, cells, losses, solve_joint, apply_rows, snr_rows  # noqa: E402

TAUS = (1.0, 0.75, 0.5, 0.25); BLK = 64
P0 = np.exp(-np.log(U1) ** 2 / (2 * SIG ** 2)); P0 = P0 / P0.sum()


class Problem:
    def __init__(self, Gc, den, ls):
        self.Gc, self.den, self.ls, self.lsm = Gc, den, ls, ls.mean((0, 1))

    def ev(self, lam):
        pi = P0 + lam; H = (self.Gc @ pi) / (self.den @ pi); R = losses(H, self.Gc, self.den, self.ls).mean((0, 1))
        return float((R + self.lsm) @ P0), R


def solve_cb(P, B, lam0, its=800):
    lam = lam0.copy(); F0, R = P.ev(lam); g = F0 + lam @ (R - B); eta = 10.0
    for it in range(its):
        new = np.maximum(0.0, lam + eta * (R - B))
        if np.array_equal(new, lam):
            break
        F2, R2 = P.ev(new); g2 = F2 + new @ (R2 - B)
        if g2 >= g - 1e-15:
            lam, F0, R, g = new, F2, R2, g2; eta = min(eta * 1.5, 1e6)
        else:
            eta *= 0.3
        if max(0.0, float((R - B).max())) <= 1e-4 * B and abs(float(lam @ (R - B))) <= 1e-6 or eta < 1e-10:
            break
    return lam, dict(B=float(B), viol_rel=float(max(0.0, (R - B).max()) / B), slack=float(lam @ (R - B)), its=it + 1,
                     F0=float(F0), R_max=float(R.max()), lam=lam.tolist(), lam_sum=float(lam.sum()))


def operator(xu, lam, f, lf, s, inv):
    out = []
    for a in range(0, len(xu), BLK):
        b = slice(a, a + BLK); Gc, den, _, sc = cells(xu[b], f, lf, s, inv[b]); pi = P0 + lam
        out.append((Gc @ pi) / (den @ pi) * sc)
    return np.concatenate(out)


def design_cb(D, cm, Q, c, tw, t_anchor, inv, taus=TAUS, g_db=G_IQ):
    """{tau: output} and the info for one input Q; taus solved in the given order with warm start."""
    f, lf, x_hat, dx_hat = setup(cm, Q, c, tw, t_anchor); s = semi.gain_sigma2(g_db); xu, _ = profiles(x_hat, dx_hat, "scalar")
    Gc, den, ls, _ = cells(xu[::2], f, lf, s, inv[::2]); _, mm = solve_joint(Gc, den, ls); R_mm = mm["value"]
    P = Problem(Gc, den, ls); _, R0 = P.ev(np.zeros(len(U1))); R_T = float(R0.max())
    info = dict(R_mm=R_mm, R_mm_gap=mm["gap"], R_T=R_T, R_T_all=R0.tolist(), runs={}); outs = {}; lam = np.zeros(len(U1))
    for tau in taus:
        lam, ri = solve_cb(P, R_mm + tau * (R_T - R_mm), lam); info["runs"][str(tau)] = ri
        outs[tau] = apply_rows(operator(xu, lam, f, lf, s, inv), D)
    return outs, info


if __name__ == "__main__" and sys.argv[1:] == ["smoke"]:
    from uqiq import TW, load
    import fdcheck as fd
    from mmr import h_of
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); arr, meta = load("ZA_c0.75_r0"); D = arr["D"]
    inv = snr_rows(D, TW, np.fft.rfftfreq(NFFT, DT))
    f, lf, xh, dx = setup(cm, Q, meta["c_eff_cum"], TW, 1.4); xu, _ = profiles(xh, dx, "scalar"); s = semi.gain_sigma2(G_IQ)
    Gc, den, _, sc = cells(xu, f, lf, s, inv); a = apply_rows(h_of(P0, Gc, den) * sc, D); b = apply_rows(operator(xu, np.zeros(17), f, lf, s, inv), D)
    print("blocked operator vs one-piece Bayes form rel %.1e" % (np.abs(a - b).max() / np.abs(a).max()))
    for m in (0.5, 1.0, 2.0):
        t0 = time.monotonic(); outs, info = design_cb(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, inv)
        print("semi m %.2f  R_mm %.4f R_T %.4f  %.0fs" % (m, info["R_mm"], info["R_T"], time.monotonic() - t0))
        for t, ri in info["runs"].items():
            nz = [(round(float(u), 2), round(l, 3)) for u, l in zip(U1, ri["lam"]) if l > 1e-6]
            print("   tau %s viol_rel %.1e slack %.1e its %d lam %s finite %s" % (t, ri["viol_rel"], ri["slack"], ri["its"], nz, bool(np.isfinite(outs[float(t)]).all())))
    D0, T, cmf, tw, fp = fd.load_fd("H1"); Df = D0 + fd.noise(D0, 14.0, 1); invf = snr_rows(Df, tw, np.fft.rfftfreq(NFFT, DT))
    outs, info = design_cb(Df, cmf, cmf["Q"], 0.5, tw, 0.0, invf)
    print("FD H1 m 0.5", {t: (round(r["viol_rel"], 6), round(r["lam_sum"], 3)) for t, r in info["runs"].items()})
