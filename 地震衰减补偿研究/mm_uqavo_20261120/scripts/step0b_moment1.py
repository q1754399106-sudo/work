"""Step 0b (POST HOC, prompted by step 0's mm0 result; descriptive, development models M1/M12 only, closed form, no noise).

(1) Why does mm0 (sigma-0 moment matching) have a smaller sensitivity than the offset-consistent sat?
    At u = 1, R_x(f) = Re[G H] (net transfer). d ln Gamma(x)/d ln u = -pi a(x) f_eff(x), f_eff(x) = int W f R_x / int W R_x.
    dGA is insensitive to u to first order iff a(x) f_eff(x) is offset-independent. Report a(x) f_eff(x)/(a(0) f_eff(0)) at
    the far trace for tik, sat, mm0.
(2) First-moment matched operator mm1 (two constraints at the model Q, deterministic):
    H = G*(1 + n1 W + n2 W f)/(|G|^2 + rho),  with  int W R_x = Gamma0  and  int W f R_x = Gamma0 f_eff(0) a(0)/a(x),
    i.e. Gamma(x;1) = Gamma(0;1) and d ln Gamma/d ln u offset-independent. Also mm1u: the same with the UQ moments (sigma) in
    place of G (E_u version of both constraints at the prior mean). Report q68, u = 1 bias, sensitivity, the minimum of
    (1 + n1 W + n2 W f) over the band (sign flip), near Gamma and far noise gain vs tik, and far-trace f_eff / near f_eff.
usage: step0b_moment1.py -> results/step0b_moment1.json"""
import sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT.parent / "uq_avo_20261118" / "scripts"))
import gate1 as g, ua  # noqa: E402
from core import ray, C  # noqa: E402
from step0_mm import op_mm, dga_u, wq, E, WT  # noqa: E402

F = g.F


def op_mm1(a_x, sig, rho, S2):
    EG, EG2 = g.moments(a_x, sig); W = S2 / S2.sum(); T = np.abs(EG) ** 2 / (EG2 + rho)
    m = lambda w: (W * w * T).sum(1)
    I0, If = m(1.0), m(F)
    A11, A12, A22 = m(W), m(W * F), m(W * F ** 2)
    G0 = I0[0]; feff0 = If[0] / I0[0]; a = np.asarray(a_x, float)
    b1 = G0 - I0; b2 = G0 * feff0 * a[0] / a - If
    det = A11 * A22 - A12 ** 2; n1 = (b1 * A22 - b2 * A12) / det; n2 = (A11 * b2 - A12 * b1) / det
    mult = 1 + n1[:, None] * W[None, :] + n2[:, None] * (W * F)[None, :]
    band = S2 > 0.05 * S2.max()
    return np.conj(EG) * mult / (EG2 + rho), float(mult[:, band].min())


def feff_ratio(H, a_x, S2):
    W = S2 / S2.sum(); R = (g.G(a_x) * H).real; fe = (W * F * R).sum(1) / (W * R).sum(1)
    return float(a_x[-1] * fe[-1] / (a_x[0] * fe[0])), float(fe[-1] / fe[0])


def main():
    out = []
    for M in ("M1", "M12"):
        Q = np.array(g.MODELS[M], float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
        for fp in g.FPS:
            S2 = g.spec(fp) ** 2
            for lam in (1e-4, 10 ** -3.5, 1e-3):
                Ht = g.op_uq(a_ev, 0.0, lam)
                base = {"tik": Ht, "sat": ua.op_ua(a0, a_ev, 0.0, lam), "mm0": np.stack([op_mm(a_ev[j], 0.0, lam, S2)[0] for j in range(3)])}
                m1 = [op_mm1(a_ev[j], 0.0, lam, S2) for j in range(3)]; base["mm1"] = np.stack([m[0] for m in m1])
                for sig in (0.447, 0.09):
                    ops = dict(base); m1u = [op_mm1(a_ev[j], sig, lam, S2) for j in range(3)]; ops["mm1u"] = np.stack([m[0] for m in m1u])
                    us = np.exp(sig * E)
                    row = dict(model=M, fp=fp, log10lam=round(float(np.log10(lam)), 2), sigma=sig,
                               mm1_minmult=[m[1] for m in m1], mm1u_minmult=[m[1] for m in m1u])
                    for name, H in ops.items():
                        pw = np.array([[dga_u(H[j], a_ev[j], s2[j], co[j], u, S2) for j in range(3)] for u in us])
                        e = 0.01; sens = [(dga_u(H[j], a_ev[j], s2[j], co[j], np.exp(e), S2) - dga_u(H[j], a_ev[j], s2[j], co[j], np.exp(-e), S2)) / (2 * e) for j in range(3)]
                        fr = [feff_ratio(H[j], a_ev[j], S2) for j in range(3)]
                        row[name] = dict(q68=float(np.median([wq(np.abs(pw[:, j]), WT, 0.68) for j in range(3)])),
                                         u1=[dga_u(H[j], a_ev[j], s2[j], co[j], 1.0, S2) for j in range(3)], sens=sens,
                                         af_ratio=[r[0] for r in fr], feff_far_over_near=[r[1] for r in fr],
                                         near=[float(g.gamma_closed(a_ev[j][:1], 1.0, H[j][:1], S2)[0] / g.gamma_closed(a_ev[j][:1], 1.0, Ht[j][:1], S2)[0]) for j in range(3)],
                                         far_noise=[float((np.abs(H[j][-1]) ** 2).sum() / (np.abs(Ht[j][-1]) ** 2).sum()) for j in range(3)])
                    out.append(row)
                    print(M, fp, "lam %.1e sig %.3f |" % (lam, sig), "  ".join(
                        "%s q68 %.3f s %+.2f af %.2f fe %.2f nz %.1f" % (k, row[k]["q68"], np.median(row[k]["sens"]), np.median(row[k]["af_ratio"]),
                                                                     np.median(row[k]["feff_far_over_near"]), np.median(row[k]["far_noise"])) for k in ops),
                          "| minmult mm1 %.2f mm1u %.2f" % (min(row["mm1_minmult"]), min(row["mm1u_minmult"])), flush=True)
    (ROOT / "results" / "step0b_moment1.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
