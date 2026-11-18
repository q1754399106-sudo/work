"""Development only (used models M1, M12): closed-form dGA of UQ-A vs Tikhonov, UQ, OC0 over u, aligned, no noise."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate1 as g, ua
from core import ray, C

U = [0.5, 0.7, 1.0, 1.4, 2.0]
for M in ["M1", "M12"]:
    Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q)
    a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T  # (3, NX) zero-offset vertical t*
    for fp in g.FPS:
        S2 = g.spec(fp) ** 2
        for lam in (10 ** -4.0, 10 ** -3.5, 10 ** -3.0):
            ops = {"tik": g.op_uq(a_ev, 0.0, lam), "uq": g.op_uq(a_ev, 0.447, lam),
                   "oc0": np.stack([g.op_oc(a_ev[j][:, None], 0.0, lam)[:, 0] for j in range(3)]),
                   "ua_nosat": ua.op_ua(a0, a_ev, 0.447, lam, sat=False), "ua": ua.op_ua(a0, a_ev, 0.447, lam),
                   "ua0": ua.op_ua(a0, a_ev, 0.0, lam)}
            line = []
            for name, H in ops.items():
                v = [[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in U]
                v = np.array(v); line.append("%s %s" % (name, " ".join("%+.2f" % x for x in np.median(v, 1))))
            print(M, fp, "lam %.0e |" % lam, " | ".join(line))
    # noise gain of ua at far offset vs zero offset (lam 10^-3.5, sat)
    lam = 10 ** -3.5; H = ua.op_ua(a0, a_ev, 0.447, lam); H0 = ua.op_ua(a0, a0, 0.447, lam)
    print(M, "noise gain far/zero per event:", np.round(np.sum(np.abs(H[:, -1]) ** 2, -1) / np.sum(np.abs(H0[:, 0]) ** 2, -1), 2))
