"""Restricted-Bayes compensation with the low-attenuation tail kept (PLAN_CBL_ZH.md):

    min_H  mean_c sum_k p_k l_ck(H)   (36 profiles, u_k = 2 * 4^(-k/16), k = 0..35, down to about 0.097)
    s.t.   R_j(H) <= B_tau             (the 17 profiles with u in [0.5, 2], i.e. mmr.U1)

p_k proportional to exp(-ln^2 u_k / 2 sigma^2) (sigma 0.447): the lognormal truncated only at u <= 2 (large u gives
negative attenuation above the semi-synthetic anchor; small u stays between x_hat and x_hat(anchor)).
B_tau = R_mm + tau (R_T - R_mm), R_mm the minimax-regret value on the 17 profiles, R_T = max_j R_j at lambda = 0 (UQWl).
Bayes form with prior p + lambda, lambda >= 0 only on the 17 cap profiles; solver cb.solve_cb (projected gradient).
Regrets do not depend on the cell normalisation, so the cap and R_mm use the 17-profile slice of the 36-profile cells.

usage: cbl.py smoke    equivalence, constraint and timing checks only (no scores)
"""
import sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mmr import np, semi, DT, NFFT, SIG, G_IQ, U1, setup, cells, losses, solve_joint, apply_rows, snr_rows  # noqa: E402
from cb import solve_cb, BLK  # noqa: E402

UL = np.r_[2 * 4.0 ** (-np.arange(35, 16, -1) / 16), U1]          # 19 lower points (ascending) + the 17-point grid
CON = np.arange(19, 36); TAUS = (1.0, 0.5, 0.25)
PL = np.exp(-np.log(UL) ** 2 / (2 * SIG ** 2)); PL = PL / PL.sum()


def prof(x_hat, dx_hat, u):
    return x_hat[:, None] + dx_hat[:, None] * (u - 1)[None, :]


class Problem:
    """ev(lam) -> (objective, regrets on the cap profiles); lam lives on the cap profiles only."""
    def __init__(self, Gc, den, ls, p, con):
        self.Gc, self.den, self.ls, self.lsm, self.p, self.con = Gc, den, ls, ls.mean((0, 1)), p, con

    def prior(self, lam):
        pi = self.p.copy(); pi[self.con] += lam; return pi

    def ev(self, lam):
        pi = self.prior(lam); H = (self.Gc @ pi) / (self.den @ pi); R = losses(H, self.Gc, self.den, self.ls).mean((0, 1))
        return float((R + self.lsm) @ self.p), R[self.con]


def operator(xu, pi, f, lf, s, inv):
    out = []
    for a in range(0, len(xu), BLK):
        b = slice(a, a + BLK); Gc, den, _, sc = cells(xu[b], f, lf, s, inv[b]); out.append((Gc @ pi) / (den @ pi) * sc)
    return np.concatenate(out)


def design_cbl(D, cm, Q, c, tw, t_anchor, inv, taus=TAUS, g_db=G_IQ, u=UL, p=PL, con=CON):
    """{tau: output} and the info for one input Q; taus solved in the given order with warm start."""
    f, lf, x_hat, dx_hat = setup(cm, Q, c, tw, t_anchor); s = semi.gain_sigma2(g_db); xu = prof(x_hat, dx_hat, u)
    Gc, den, ls, _ = cells(xu[::2], f, lf, s, inv[::2])
    _, mm = solve_joint(Gc[..., con], den[..., con], ls[..., con]); R_mm = mm["value"]
    P = Problem(Gc, den, ls, p, con); _, R0 = P.ev(np.zeros(len(con))); R_T = float(R0.max())
    info = dict(R_mm=R_mm, R_mm_gap=mm["gap"], R_T=R_T, R_T_all=R0.tolist(), runs={}); outs = {}; lam = np.zeros(len(con))
    for tau in taus:
        lam, ri = solve_cb(P, R_mm + tau * (R_T - R_mm), lam); info["runs"][str(tau)] = ri
        outs[tau] = apply_rows(operator(xu, P.prior(lam), f, lf, s, inv), D)
    return outs, info


if __name__ == "__main__" and sys.argv[1:] == ["smoke"]:
    from uqiq import TW, load
    import fdcheck as fd
    import cb
    from mmr import GH, h_of
    print("grid", UL.size, "u_min %.4f" % UL[0], "upper 17 == U1:", bool(np.array_equal(UL[CON], U1)),
          "prior mass below 0.5 %.3f" % PL[UL < 0.5 - 1e-12].sum())
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); arr, meta = load("ZA_c0.75_r0"); D = arr["D"]
    inv = snr_rows(D, TW, np.fft.rfftfreq(NFFT, DT)); c = meta["c_eff_cum"]
    # (a) restricted to the 17-point grid this code is cb.py
    t0 = time.monotonic(); o17, i17 = design_cbl(D, cm, Q, c, TW, 1.4, inv, taus=(1.0, 0.5), u=U1, p=cb.P0, con=np.arange(17))
    oc, ic = cb.design_cb(D, cm, Q, c, TW, 1.4, inv, taus=(1.0, 0.5))
    print("17-grid vs cb.py: UQWt rel %.1e, tau 0.5 rel %.1e, R_mm %.2e, R_T %.2e  (%.0fs)" % (
        np.abs(o17[1.0] - oc[1.0]).max() / np.abs(oc[1.0]).max(), np.abs(o17[0.5] - oc[0.5]).max() / np.abs(oc[0.5]).max(),
        abs(i17["R_mm"] - ic["R_mm"]), abs(i17["R_T"] - ic["R_T"]), time.monotonic() - t0))
    # (b) feasibility, lambda support, finiteness, timing; (c) UQWl vs UQW operator difference (descriptive)
    def uqw_out(D, cm, Q, c, tw, ta, inv):
        f, lf, xh, dx = setup(cm, Q, c, tw, ta); e, w = GH; xs = xh[:, None] + dx[:, None] * (np.exp(SIG * e) - 1)[None, :]
        Gc, den, _, sc = cells(xs, f, lf, semi.gain_sigma2(G_IQ), inv); return apply_rows(h_of(w / w.sum(), Gc, den) * sc, D)
    for m in (0.5, 1.0, 2.0):
        t0 = time.monotonic(); outs, info = design_cbl(D, cm, Q, m * c, TW, 1.4, inv); dt = time.monotonic() - t0
        ref = uqw_out(D, cm, Q, m * c, TW, 1.4, inv)
        print("semi m %.2f  R_mm %.4f R_T %.4f  %.0fs  UQWl vs UQW rel %.3f" % (m, info["R_mm"], info["R_T"], dt,
              np.linalg.norm(outs[1.0] - ref) / np.linalg.norm(ref)))
        for t, ri in info["runs"].items():
            nz = [(round(float(uu), 2), round(l, 3)) for uu, l in zip(UL[CON], ri["lam"]) if l > 1e-6]
            print("   tau %s viol_rel %.1e slack %.1e its %d lam %s finite %s" % (t, ri["viol_rel"], ri["slack"], ri["its"], nz,
                                                                            bool(np.isfinite(outs[float(t)]).all())))
    D0, T, cmf, tw, fp = fd.load_fd("H1"); Df = D0 + fd.noise(D0, 14.0, 1); invf = snr_rows(Df, tw, np.fft.rfftfreq(NFFT, DT))
    t0 = time.monotonic(); outs, info = design_cbl(Df, cmf, cmf["Q"], 0.5, tw, 0.0, invf); dt = time.monotonic() - t0
    ref = uqw_out(Df, cmf, cmf["Q"], 0.5, tw, 0.0, invf)
    print("FD H1 m 0.5  %.0fs  UQWl vs UQW rel %.3f" % (dt, np.linalg.norm(outs[1.0] - ref) / np.linalg.norm(ref)),
          {t: (round(r["viol_rel"], 6), round(r["lam_sum"], 3), bool(np.isfinite(outs[float(t)]).all())) for t, r in info["runs"].items()})
