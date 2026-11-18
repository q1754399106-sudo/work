"""Development only (M1, M12): prior-weighted RMS of dGA vs prior width sig, when the operator is built at the model Q
(u-median 1) and the truth is u ~ lognormal(0, sig). Closed form, aligned, no noise. tik: one-step Tikhonov;
sat~: sigma-0 offset-consistent (~ saturated damping); ud: zero-offset UQ with the same sig x deterministic offset correction."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate1 as g, ua
from dev_ud import ud
from core import ray, C

for M in ["M1", "M12"]:
    Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
    for fp in g.FPS:
        S2 = g.spec(fp) ** 2; lam = 10 ** -3.5
        for sig in (0.0, 0.05, 0.1, 0.2, 0.3, 0.447):
            un = np.exp(sig * g.GH_E); res = []
            for name, H in (("tik", g.op_uq(a_ev, 0.0, lam)), ("sat~", ua.op_ua(a0, a_ev, 0.0, lam)), ("ud", ud(a0, a_ev, sig, lam, 1.0))):
                pw = np.array([[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in un])
                res.append("%s %.3f" % (name, np.median(np.sqrt(g.GH_W @ pw ** 2))))
            print(M, fp, "sig %.2f:" % sig, "  ".join(res))
        # sensitivity d dGA / d ln u at u = 1 for sat~
        H = ua.op_ua(a0, a_ev, 0.0, lam); e = 0.05
        d = [(g.dga(g.gamma_closed(a_ev[j], np.exp(e), H[j], S2), s2[j], co[j]) - g.dga(g.gamma_closed(a_ev[j], np.exp(-e), H[j], S2), s2[j], co[j])) / (2 * e) for j in range(3)]
        print(M, fp, "sat~ d dGA / d ln u per event:", np.round(d, 2))
