"""Step 0d (POST HOC, development models M1/M12, closed form, no noise; descriptive).

Claim to test (THEORY_ZH.md section 5): an operator can be unbiased at the model Q for EVERY local reflectivity spectrum only if
its net transfer R_x(f) = Re[G H] is the same on all offsets (offset-consistent). mm0 / mm1 / fs lower the Q-level
sensitivity precisely by making far-trace R_x differ (lower effective frequency), so they must become biased at u = 1 as soon
as the scoring weight differs from the design weight. The scoring weight of a per-event amplitude read by correlation with the
true (wavelet-convolved) event is W_r = S^2 bp^2 |R(f)|^2, R(f) = local reflectivity spectrum.
Design weight for all operators: Wd = S^2 bp^2 (isolated event, exact wavelet). Scoring weights:
  iso  |R|^2 = 1;  blue05, blue10  |R|^2 = (f/fp)^0.5, (f/fp)^1;  thin8, thin16, thin24  |R|^2 = 4 sin^2(pi f tau), tau 8/16/24 ms
  (the step-2b thin scene is tau = 16 ms).
Operators: tik, sat, mm0, mm1, fs (step 0c) and fsc (fs with R_x clipped to <= max R_0, then rescaled).
Reported per weight: |dGA| at u = 1 (median over events), prior 68th percentile of |dGA| (sigma 0.09, 0.447).
usage: step0d_reflspec.py -> results/step0d_reflspec.json"""
import sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT.parent / "uq_avo_20261118" / "scripts"))
import gate1 as g, ua  # noqa: E402
from core import ray, C  # noqa: E402
from step0_mm import op_mm, dga_u, wq, E, WT  # noqa: E402
from step0b_moment1 import op_mm1  # noqa: E402
from step0c_fscale import op_fs, bandw  # noqa: E402

F = g.F


def op_fsc(a_x, rho, Wd):
    a = np.asarray(a_x, float); G0 = g.G(a[0]); R0 = np.abs(G0) ** 2 / (np.abs(G0) ** 2 + rho); H = op_fs(a_x, rho, Wd)
    for ix in range(len(a)):
        Gx = g.G(a[ix]); Rx = np.minimum((Gx * H[ix]).real, R0.max()); Rx *= (Wd * R0).sum() / max((Wd * Rx).sum(), 1e-300)
        H[ix] = np.where(np.abs(Gx) > 1e-300, Rx * np.conj(Gx) / np.maximum(np.abs(Gx) ** 2, 1e-300), 0.0)
    return H


def main():
    out = []
    for M in ("M1", "M12"):
        Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
        for fp in g.FPS:
            Wd = g.spec(fp) ** 2 * bandw(fp)
            weights = {"iso": Wd, "blue05": Wd * (F / fp) ** 0.5, "blue10": Wd * (F / fp),
                       **{"thin%d" % int(tau * 1000): Wd * 4 * np.sin(np.pi * F * tau) ** 2 for tau in (0.008, 0.016, 0.024)}}
            for lam in (1e-4, 10 ** -3.5):
                ops = {"tik": g.op_uq(a_ev, 0.0, lam), "sat": ua.op_ua(a0, a_ev, 0.0, lam),
                       "mm0": np.stack([op_mm(a_ev[j], 0.0, lam, Wd)[0] for j in range(3)]),
                       "mm1": np.stack([op_mm1(a_ev[j], 0.0, lam, Wd)[0] for j in range(3)]),
                       "fs": np.stack([op_fs(a_ev[j], lam, Wd) for j in range(3)]), "fsc": np.stack([op_fsc(a_ev[j], lam, Wd) for j in range(3)])}
                row = dict(model=M, fp=fp, log10lam=round(float(np.log10(lam)), 2))
                for wn, Ws in weights.items():
                    row[wn] = {}
                    for name, H in ops.items():
                        d = lambda j, u: dga_u(H[j], a_ev[j], s2[j], co[j], u, Ws)
                        r = dict(u1=float(np.median([abs(d(j, 1.0)) for j in range(3)])), u1_signed=[d(j, 1.0) for j in range(3)])
                        for sig in (0.09, 0.447):
                            pw = np.array([[d(j, u) for j in range(3)] for u in np.exp(sig * E)])
                            r["q68_%g" % sig] = float(np.median([wq(np.abs(pw[:, j]), WT, 0.68) for j in range(3)]))
                        row[wn][name] = r
                    print(M, fp, "lam %.1e %-6s |" % (lam, wn), "  ".join("%s u1 %.3f q09 %.3f q45 %.3f" % (k, v["u1"], v["q68_0.09"], v["q68_0.447"]) for k, v in row[wn].items()), flush=True)
                out.append(row)
    (ROOT / "results" / "step0d_reflspec.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
