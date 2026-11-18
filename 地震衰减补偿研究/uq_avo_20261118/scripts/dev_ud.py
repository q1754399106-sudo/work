"""Development only (M1, M12): decoupled UQ  H = H0_UQ(a0) * Rdet,  Rdet = 1/G(ubar * ax) (saturated, gamma 6),
against Tikhonov, UQ, the sigma-0 offset-consistent operator (~ saturated damping) and UQ-A. Closed form, aligned, no noise.
Prior-weighted RMS and mean of dGA over u ~ lognormal(0, 0.447) (41-pt GH), and dGA on the u grid."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate1 as g, ua
from core import ray, C

SIG = 0.447; U = [0.5, 0.7, 1.0, 1.4, 2.0]; UQN = np.exp(SIG * g.GH_E); W = g.GH_W


def ud(a0, a, sig, rho, ubar):
    E0, E02 = g.moments(a0, sig); H0 = np.conj(E0) / (E02 + rho); R = 1.0 / g.G(ubar * (a - a0))
    R, _ = ua.saturate(H0, R); return H0 * R


for M in ["M1", "M12"]:
    Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
    for fp in g.FPS:
        S2 = g.spec(fp) ** 2
        for lam in (10 ** -4.0, 10 ** -3.5):
            ops = {"tik": g.op_uq(a_ev, 0.0, lam), "uq": g.op_uq(a_ev, SIG, lam), "sat~": ua.op_ua(a0, a_ev, 0.0, lam),
                   "uqa": ua.op_ua(a0, a_ev, SIG, lam), "ud1": ud(a0, a_ev, SIG, lam, 1.0), "udm": ud(a0, a_ev, SIG, lam, np.exp(SIG ** 2 / 2))}
            out = []
            for name, H in ops.items():
                grid = np.median([[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in U], 1)
                pw = np.array([[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in UQN])  # (q, 3)
                rms = np.sqrt(W @ pw ** 2); mean = W @ pw
                out.append("%s rms %.2f mean %+.2f [%s]" % (name, np.median(rms), np.median(mean), " ".join("%+.2f" % x for x in grid)))
            print(M, fp, "lam %.0e" % lam); print("   " + "\n   ".join(out))
    # ubar that minimises prior RMS (scan), lam 10^-3.5, 40 and 60 Hz
    lam = 10 ** -3.5
    for fp in g.FPS:
        S2 = g.spec(fp) ** 2; best = []
        for ub in np.linspace(0.8, 1.4, 13):
            H = ud(a0, a_ev, SIG, lam, ub); pw = np.array([[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in UQN])
            best.append((np.median(np.sqrt(W @ pw ** 2)), ub))
        print(M, fp, "rms vs ubar:", " ".join("%.2f@%.2f" % b for b in best))
