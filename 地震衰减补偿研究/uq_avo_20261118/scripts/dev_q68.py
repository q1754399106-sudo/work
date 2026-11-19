"""Robustness check of DEV_GATE2_ZH.md: the prior-weighted RMS of dGA is tail-dominated at sigma >= 0.35 (G/A is a ratio that
diverges when the compensated intercept nears zero). Recompute with the prior 68th percentile of |dGA| (truncated +-3 sigma,
121 points), for the operators of dev_ud.py, M1 and M12, lambda 10^-3.5."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate1 as g, ua
from dev_ud import ud
from core import ray, C


def wq(x, w, q):
    i = np.argsort(x); cw = np.cumsum(w[i]); return float(x[i][np.searchsorted(cw, q * cw[-1])])


E = np.linspace(-3, 3, 121); Wt = np.exp(-E ** 2 / 2); Wt /= Wt.sum(); out = []; lam = 10 ** -3.5
for M in ("M1", "M12"):
    Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
    for fp in g.FPS:
        S2 = g.spec(fp) ** 2
        for sig in (0.1, 0.2, 0.3, 0.447):
            ops = {"tik": g.op_uq(a_ev, 0.0, lam), "uq": g.op_uq(a_ev, sig, lam), "sat~": ua.op_ua(a0, a_ev, 0.0, lam),
                   "uqa": ua.op_ua(a0, a_ev, sig, lam), "ud1": ud(a0, a_ev, sig, lam, 1.0)}
            row = dict(model=M, fp=fp, sigma=sig)
            for name, H in ops.items():
                pw = np.array([[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in np.exp(sig * E)])
                row[name] = float(np.median([wq(np.abs(pw[:, j]), Wt, 0.68) for j in range(3)]))
            out.append(row); print(M, fp, sig, "  ".join("%s %.3f" % (k, row[k]) for k in ops), flush=True)
(Path(__file__).resolve().parents[1] / "results" / "dev_q68.json").write_text(json.dumps(out, indent=1))
