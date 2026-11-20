"""Step 0 of mm_uqavo_20261120 (THEORY_ZH.md): can the moment-matched UQ-AVO operator beat Tikhonov or spectral-ratio calibration?

Closed-form Gamma only (validated in uq_avo_20261118 gate 1, P1: slope 0.974, R^2 0.999), aligned gather, no noise,
development models M1 and M12 only (already used; no confirmation data consumed). Descriptive, no gate.

Operators per event j and offset x (row at the event time, as in uq_avo dev_q68.py):
  tik    one-step Tikhonov                      G*/(|G|^2+rho)
  uq     UQ-W gamma 0                           E[G]*/(E|G|^2+rho)
  oc     offset-consistent UQ (gate 1 OC)       equalise E-transfer to the far-offset minimum, per frequency
  sat    sigma-0 offset-consistent (~ saturated damping; ua.op_ua sig 0)
  uqa    UQ-A, prior-expected response ratio (ua.op_ua sig)
  mm     moment-matched UQ (user's proposal)    E[G]*(1+nu W)/(E|G|^2+rho), W = S^2/sum S^2,
                                                nu(x) = (Gamma0 - I1(x))/I2(x), Gamma0 = I1(0), I_k = sum W^k |E G|^2/(E|G|^2+rho)
  mm0    same with sigma = 0 (constraint on Gamma at u = 1 instead of E_u Gamma)
  mmsat  saturated-damping kernel + moment matching:  H_sat (1+nu W), nu s.t. E_u Gamma(x) = E_u Gamma(0)
Metrics: dGA on u grid; prior mean and 68th percentile of |dGA| (u ~ lognormal(0, sigma_true), +-3 sigma, 121 pts);
first-order sensitivity s = d dGA / d ln u at u = 1; near-trace Gamma(0, u=1) and far-trace noise gain sum|H|^2, both / tik.
Scenarios: (A) sigma_op = sigma_true = 0.447 (no calibration); (B) sigma_op = sigma_true = 0.09 (after offset spectral-ratio
calibration, ln(c_hat/c) RMS 0.069-0.096 in offset_q_avo gates A2/A3).
usage: step0_mm.py -> results/step0_mm.json"""
import os, sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[0]
sys.path.insert(0, str(ROOT.parent / "uq_avo_20261118" / "scripts"))
import gate1 as g, ua  # noqa: E402
from core import ray, C  # noqa: E402

E = np.linspace(-3, 3, 121); WT = np.exp(-E ** 2 / 2); WT /= WT.sum()


def wq(x, w, q):
    i = np.argsort(x); cw = np.cumsum(w[i]); return float(x[i][np.searchsorted(cw, q * cw[-1])])


def op_mm(a_x, sig, rho, S2):
    """a_x (NX,) path t* of one event over offsets -> (H (NX, nf), nu (NX,))."""
    EG, EG2 = g.moments(a_x, sig); W = S2 / S2.sum(); T = np.abs(EG) ** 2 / (EG2 + rho)
    I1 = (W * T).sum(1); I2 = (W ** 2 * T).sum(1); nu = (I1[0] - I1) / I2
    return np.conj(EG) * (1 + nu[:, None] * W[None, :]) / (EG2 + rho), nu


def op_mmsat(a0, a_x, sig, rho, S2):
    """saturated-damping kernel (sigma 0) with a band-integrated moment-matching correction on E_u Gamma."""
    Hs = ua.op_ua(a0, a_x, 0.0, rho); EG, _ = g.moments(a_x, sig); W = S2 / S2.sum()
    R = (EG * Hs).real; J1 = (W * R).sum(1); J2 = (W ** 2 * R).sum(1); nu = (J1[0] - J1) / J2
    return Hs * (1 + nu[:, None] * W[None, :]), nu


def dga_u(Hj, a_j, s2j, coj, u, S2):
    return g.dga(g.gamma_closed(a_j, u, Hj, S2), s2j, coj)


def main():
    out = []
    for M in ("M1", "M12"):
        Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
        for fp in g.FPS:
            S2 = g.spec(fp) ** 2
            for lam in (1e-4, 10 ** -3.5, 1e-3):
                for sig in (0.447, 0.09):
                    ops = {"tik": g.op_uq(a_ev, 0.0, lam), "uq": g.op_uq(a_ev, sig, lam),
                           "oc": np.stack([g.op_oc(a_ev[j][:, None], sig, lam)[:, 0] for j in range(3)]),
                           "sat": ua.op_ua(a0, a_ev, 0.0, lam), "uqa": ua.op_ua(a0, a_ev, sig, lam)}
                    mm = [op_mm(a_ev[j], sig, lam, S2) for j in range(3)]; ops["mm"] = np.stack([m[0] for m in mm])
                    mm0 = [op_mm(a_ev[j], 0.0, lam, S2) for j in range(3)]; ops["mm0"] = np.stack([m[0] for m in mm0])
                    ms = [op_mmsat(a0[j], a_ev[j], sig, lam, S2) for j in range(3)]; ops["mmsat"] = np.stack([m[0] for m in ms])
                    us = np.exp(sig * E); row = dict(model=M, fp=fp, log10lam=round(float(np.log10(lam)), 2), sigma=sig,
                                                       nu_far=[float(m[1][-1]) for m in mm], nu_far_sat=[float(m[1][-1]) for m in ms])
                    Ht = ops["tik"]
                    for name, H in ops.items():
                        pw = np.array([[dga_u(H[j], a_ev[j], s2[j], co[j], u, S2) for j in range(3)] for u in us])
                        grid = [[dga_u(H[j], a_ev[j], s2[j], co[j], u, S2) for j in range(3)] for u in (0.5, 0.7, 1.0, 1.4, 2.0)]
                        e = 0.01; sens = [(dga_u(H[j], a_ev[j], s2[j], co[j], np.exp(e), S2) - dga_u(H[j], a_ev[j], s2[j], co[j], np.exp(-e), S2)) / (2 * e) for j in range(3)]
                        near = [float(g.gamma_closed(a_ev[j][:1], 1.0, H[j][:1], S2)[0] / g.gamma_closed(a_ev[j][:1], 1.0, Ht[j][:1], S2)[0]) for j in range(3)]
                        noise = [float((np.abs(H[j][-1]) ** 2).sum() / (np.abs(Ht[j][-1]) ** 2).sum()) for j in range(3)]
                        row[name] = dict(q68=float(np.median([wq(np.abs(pw[:, j]), WT, 0.68) for j in range(3)])),
                                         q68_ev=[wq(np.abs(pw[:, j]), WT, 0.68) for j in range(3)],
                                         mean=float(np.median(WT @ pw)), u1=float(np.median(np.abs(grid[2]))), u1_signed=grid[2],
                                         grid_med=[float(np.median(r)) for r in grid], sens=sens, near=near, far_noise=noise)
                    out.append(row)
                    print(M, fp, "lam %.1e sig %.3f |" % (lam, sig), "  ".join("%s q68 %.3f u1 %+.3f s %+.2f" % (k, row[k]["q68"], np.median(row[k]["u1_signed"]), np.median(row[k]["sens"])) for k in ops), flush=True)
    (ROOT / "results" / "step0_mm.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
