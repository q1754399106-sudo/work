"""Data for paper figure 3 (same settings as dev_sigma.py, finer sigma grid, saved): prior-weighted RMS of per-event dGA vs
sigma_lnQ for one-step Tikhonov and the sigma-0 offset-consistent operator, M1 and M12, 60/40 Hz, lambda 10^-3.5, closed form.
Prior quadrature: e on 121 equally spaced points in [-3, 3], normal weights renormalised (truncated at 3 sigma). The 41-point
Gauss-Hermite rule of dev_sigma.py puts nodes at |e| up to ~8.6 (u up to e^{+-4} at sigma 0.45), where dGA is extreme, so its
RMS is tail-dominated and jumps between sigma values (seen in the first draft of figure 3); both are saved.
The RMS of a ratio metric is tail-dominated at sigma >= 0.35 even with truncation (the compensated intercept nears zero for
extreme u), so the figure uses the prior 68th percentile of |dGA| (q68), which equals k sigma in the linear regime."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate1 as g, ua
from core import ray, C

out = []; lam = 10 ** -3.5
for M in ("M1", "M12"):
    Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
    ops = {"tikhonov": g.op_uq(a_ev, 0.0, lam), "offset_consistent": ua.op_ua(a0, a_ev, 0.0, lam)}
    for fp in g.FPS:
        S2 = g.spec(fp) ** 2
        for sig in np.round(np.arange(0, 0.501, 0.05), 3):
            E = np.linspace(-3, 3, 121); Wt = np.exp(-E ** 2 / 2); Wt /= Wt.sum()
            for name, H in ops.items():
                pw = np.array([[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in np.exp(sig * E)])
                pg = np.array([[g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j]) for j in range(3)] for u in np.exp(sig * g.GH_E)])
                q68 = [float(np.abs(pw[:, j])[np.argsort(np.abs(pw[:, j]))][np.searchsorted(np.cumsum(Wt[np.argsort(np.abs(pw[:, j]))]), 0.68)]) for j in range(3)]
                out.append(dict(model=M, fp=fp, sigma=float(sig), op=name, q68_median=float(np.median(q68)), rms_median=float(np.median(np.sqrt(Wt @ pw ** 2))),
                                rms_events=np.sqrt(Wt @ pw ** 2).tolist(), rms_median_gh41=float(np.median(np.sqrt(g.GH_W @ pg ** 2)))))
(Path(__file__).resolve().parents[1] / "results" / "fig_sigma.json").write_text(json.dumps(out, indent=1)); print(len(out), "rows")
