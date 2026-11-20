"""Step 0c (POST HOC, development models M1/M12 only, closed form, no noise; descriptive).

Step 0b's mm1 removes the first-order Q-level sensitivity but its multiplier (1 + n1 W + n2 W f) goes negative in band
(polarity flips). Exact, sign-safe version of the same idea:
  Gamma(x;u)/Gamma(x;1) = int p_x(f) exp[-pi f (u-1) a(x)] df,  p_x = W R_x / int W R_x   (amplitude part of G; R_x = Re[G H] at u=1)
  is offset-independent for every u  iff  p_x(f) = r p_0(r f),  r = a(x)/a(0)   ("f t* invariance").
  => R_x(f) = Gamma0 r p_0(r f) / W(f) = r W(r f) R_0(r f) / W(f),  H_x = R_x / G(a(x)) (R_x >= 0, a low-pass for r > 1).
  Base near-trace transfer R_0 = Tikhonov |G0|^2/(|G0|^2 + rho). The dispersion phase is not scale-invariant (ln f term) and
  is left as is.
Weights: W = S^2 x band-pass^2 (10 Hz - fp, as the step-2b G/A metric; theory_k bandw) for both design and scoring.
Operators: tik, sat (sigma-0 offset-consistent), mm0, mm1 (step 0b), fs (frequency scaled).
Robustness (fs only): designed with a wrong wavelet (fp x 0.8, x 1.2) or without the band-pass, scored with the true W.
Metrics as step 0: prior 68th percentile of |dGA| (sigma 0.447 and 0.09), u = 1 bias, sensitivity, u grid, near-trace Gamma
and far-trace noise gain sum|H|^2 relative to tik, far/near effective frequency, far-trace sum|H|^2 relative to sat.
usage: step0c_fscale.py -> results/step0c_fscale.json"""
import sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT.parent / "uq_avo_20261118" / "scripts"))
import gate1 as g, ua  # noqa: E402
from core import ray, C  # noqa: E402
from step0_mm import op_mm, dga_u, wq, E, WT  # noqa: E402
from step0b_moment1 import op_mm1  # noqa: E402

F = g.F


def bandw(fp):
    b0, b1 = 10.0, fp; return (np.clip((F - b0 + 5) / 5, 0, 1) * np.clip((b1 + 8 - F) / 8, 0, 1)) ** 2


def op_fs(a_x, rho, Wd):
    """a_x (NX,) model path t*; Wd design weight (nf,). Returns H (NX, nf)."""
    a = np.asarray(a_x, float); G0 = g.G(a[0]); R0 = np.abs(G0) ** 2 / (np.abs(G0) ** 2 + rho); P0 = Wd * R0
    H = np.zeros((len(a), len(F)), complex); ok = Wd > 1e-12 * Wd.max()
    for ix in range(len(a)):
        r = a[ix] / a[0]; Pr = np.interp(r * F, F, P0, right=0.0) * r
        Rx = np.where(ok, Pr / np.where(ok, Wd, 1.0), 0.0)
        Rx *= (P0.sum()) / max((Wd * Rx).sum(), 1e-300)  # keep Gamma(x;1) = Gamma(0;1) under the design weight
        Gx = g.G(a[ix]); H[ix] = np.where(np.abs(Gx) > 1e-300, Rx * np.conj(Gx) / np.maximum(np.abs(Gx) ** 2, 1e-300), 0.0)
    return H


def main():
    out = []
    for M in ("M1", "M12"):
        Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
        for fp in g.FPS:
            Wt = g.spec(fp) ** 2 * bandw(fp)
            for lam in (1e-4, 10 ** -3.5, 1e-3):
                Ht = g.op_uq(a_ev, 0.0, lam); Hs = ua.op_ua(a0, a_ev, 0.0, lam)
                ops = {"tik": Ht, "sat": Hs, "mm0": np.stack([op_mm(a_ev[j], 0.0, lam, Wt)[0] for j in range(3)]),
                       "mm1": np.stack([op_mm1(a_ev[j], 0.0, lam, Wt)[0] for j in range(3)]),
                       "fs": np.stack([op_fs(a_ev[j], lam, Wt) for j in range(3)]),
                       "fs_w08": np.stack([op_fs(a_ev[j], lam, g.spec(0.8 * fp) ** 2 * bandw(fp)) for j in range(3)]),
                       "fs_w12": np.stack([op_fs(a_ev[j], lam, g.spec(1.2 * fp) ** 2 * bandw(fp)) for j in range(3)]),
                       "fs_nobp": np.stack([op_fs(a_ev[j], lam, g.spec(fp) ** 2) for j in range(3)])}
                for sig in (0.447, 0.09):
                    us = np.exp(sig * E); row = dict(model=M, fp=fp, log10lam=round(float(np.log10(lam)), 2), sigma=sig)
                    for name, H in ops.items():
                        d = lambda j, u: dga_u(H[j], a_ev[j], s2[j], co[j], u, Wt)
                        pw = np.array([[d(j, u) for j in range(3)] for u in us]); e = 0.01
                        R = [(g.G(a_ev[j]) * H[j]).real for j in range(3)]
                        fe = [float(((Wt * F * R[j][-1]).sum() / (Wt * R[j][-1]).sum()) / ((Wt * F * R[j][0]).sum() / (Wt * R[j][0]).sum())) for j in range(3)]
                        row[name] = dict(q68=float(np.median([wq(np.abs(pw[:, j]), WT, 0.68) for j in range(3)])),
                                         u1=[d(j, 1.0) for j in range(3)], sens=[(d(j, np.exp(e)) - d(j, np.exp(-e))) / (2 * e) for j in range(3)],
                                         grid=[[d(j, u) for j in range(3)] for u in (0.5, 0.7, 1.0, 1.4, 2.0)], feff_far_over_near=fe,
                                         near=[float(g.gamma_closed(a_ev[j][:1], 1.0, H[j][:1], Wt)[0] / g.gamma_closed(a_ev[j][:1], 1.0, Ht[j][:1], Wt)[0]) for j in range(3)],
                                         far_noise_vs_tik=[float((np.abs(H[j][-1]) ** 2).sum() / (np.abs(Ht[j][-1]) ** 2).sum()) for j in range(3)],
                                         far_noise_vs_sat=[float((np.abs(H[j][-1]) ** 2).sum() / (np.abs(Hs[j][-1]) ** 2).sum()) for j in range(3)],
                                         min_R_inband=float(min(R[j][:, Wt > 0.05 * Wt.max()].min() for j in range(3))))
                    out.append(row)
                    print(M, fp, "lam %.1e sig %.3f |" % (lam, sig), "  ".join(
                        "%s %.3f/%+.2f/fe%.2f/nzs%.2f" % (k, row[k]["q68"], np.median(row[k]["sens"]), np.median(row[k]["feff_far_over_near"]),
                                                       np.median(row[k]["far_noise_vs_sat"])) for k in ops), flush=True)
    (ROOT / "results" / "step0c_fscale.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
