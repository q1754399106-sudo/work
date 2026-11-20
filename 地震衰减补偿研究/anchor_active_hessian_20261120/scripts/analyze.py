"""Scores of PLAN_ZH.md section 3 from results/*.json -> analysis/hessian.json and a printed table.

Per profile: medians over datasets of every alpha, errors against K_ref (11-19 plan table), verdicts M1-M4, the exact
three-block split of J_vv and J_v va (theta-theta, theta-mu, mu-mu), coupling absorption and curvature location by
0.1 s block. a0/a(t_med) is listed as post hoc (not judged).
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]; RES = HERE / "results"; OUT = HERE / "analysis"
REF = {"base": (0.60, 0.664), "S": (1.00, 0.856), "D": (4.65, 0.598), "D2": (2.28, 0.624), "M": (0.24, 0.614),
       "R1": (1.00, 0.649), "R2": (0.42, 0.689), "R3": (0.28, 0.608), "R4": (0.38, 0.737), "X1": (1.64, 0.638),
       "X2": (3.96, 0.611), "X3": (1.39, 0.615)}
ORDER = list(REF)
R_LE_1 = [p for p in ORDER if REF[p][0] <= 1.0]


def blocks3(r):
    """Exact theta-theta / theta-mu / mu-mu split (DERIVATION_ZH.md section 5): curvature C and numerator N = -J_v va."""
    f, mp, md = r["F"]["fix"], r["mu_partials"], r["mu_derivs"]
    C = dict(tt=f["Jvv"], tm=2 * mp["J_vmu"] * md["mu_v"], mm=mp["J_mumu"] * md["mu_v"] ** 2 + mp["J_mu"] * md["mu_vv"])
    N = dict(tt=-f["Jvx"], tm=-(mp["J_vmu"] * md["mu_va"] + mp["J_vamu"] * md["mu_v"]),
             mm=-(mp["J_mumu"] * md["mu_v"] * md["mu_va"] + mp["J_mu"] * md["mu_vva"]))
    tot_C, tot_N = sum(C.values()), sum(N.values())
    return dict(C=C, N=N, alpha_block={k: N[k] / C[k] for k in C}, share_C={k: C[k] / tot_C for k in C},
                alpha_sum=tot_N / tot_C, closure=tot_C / r["F"]["pipe"]["Jvv"] - 1)


def per_dataset(r):
    al = dict(r["alpha"]); a0 = r["a0"]
    al.update(R=r["R"]["0.05"]["alpha"], F=r["F"]["pipe"]["alpha"], F_fix=r["F"]["fix"]["alpha"],
              F_theta=r["F"]["fix_theta"]["alpha"], F_amp=r["F"]["fix_amp"]["alpha"],
              F_q25=r["F"]["q25"]["alpha"], F_q75=r["F"]["q75"]["alpha"], med=a0 / r["a_at"]["t_med"],
              med_q25=a0 / r["a_at"]["q25"], med_q75=a0 / r["a_at"]["q75"])
    return al


def depth_tables(rs):
    """0.1 s blocks, summed over the datasets of a profile: coupling absorption 1 - diag H / diag Hd, and the share of
    the theta-theta curvature a'Ha (diagonal blocks of a o H o a) against the energy weight share."""
    t = np.asarray(rs[0]["t"]); edges = np.asarray(rs[0]["blocks01"]["edges"]); nb = len(edges) - 1
    seg = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, nb - 1)
    dH = np.zeros(nb); dHd = np.zeros(nb); den = np.zeros(nb); denD = np.zeros(nb); en = np.zeros(nb); dens = np.zeros(nb)
    for r in rs:
        a = np.asarray(r["a"]); dp = r["depth"]
        for b in range(nb):
            m = seg == b
            dH[b] += np.sum(np.asarray(dp["diagH"])[m]); dHd[b] += np.sum(np.asarray(dp["diagHd"])[m])
            den[b] += np.sum((a * np.asarray(dp["Ha"]))[m]) / r["curv"]["aHa"]
            denD[b] += np.sum((a * np.asarray(dp["Hda"]))[m]) / r["curv"]["aHda"]
            e = np.asarray(dp["e"]); en[b] += np.sum((e * a ** 2)[m]) / np.sum(e * a ** 2)
            dens[b] += np.mean(np.asarray(dp["n"])[m]) if m.any() else 0.0
    k = len(rs); keep = np.array([np.any(seg == b) for b in range(nb)])      # drop empty 0.1 s blocks
    dH, dHd, den, denD, en, dens = (x[keep] for x in (dH, dHd, den, denD, en, dens))
    mids = ((np.asarray(rs[0]["blocks01"]["edges"])[:-1] + np.asarray(rs[0]["blocks01"]["edges"])[1:]) / 2)[keep]
    return dict(t_mid=mids.round(2).tolist(), absorbed=(1 - dH / dHd).round(3).tolist(),
                share_aHa=(den / k).round(3).tolist(), share_aHda=(denD / k).round(3).tolist(),
                share_energy_a2=(en / k).round(3).tolist(), active_density=(dens / k).round(3).tolist())


def main():
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(RES.glob("*.json"))]
    groups = {}
    for r in recs:
        groups.setdefault(r["profile"], []).append(r)
    out = dict(profiles={}, verdicts={}, n_datasets=len(recs))
    keys = None; dF = []; inter = []
    for p in ORDER:
        rs = groups.get(p, [])
        if not rs:
            continue
        ds = [per_dataset(r) for r in rs]; keys = list(ds[0])
        bl = [blocks3(r) for r in rs]
        rr, K = REF[p]
        row = dict(r=rr, K_ref=K, n=len(rs), all_converged=all(r["all_converged"] for r in rs),
                   alpha={k: float(np.median([d[k] for d in ds])) for k in keys},
                   alpha_range={k: [float(min(d[k] for d in ds)), float(max(d[k] for d in ds))] for k in keys},
                   blocks=dict(share_C={k: float(np.median([b["share_C"][k] for b in bl])) for k in ("tt", "tm", "mm")},
                               alpha_block={k: float(np.median([b["alpha_block"][k] for b in bl])) for k in ("tt", "tm", "mm")},
                               closure_max=float(max(abs(b["closure"]) for b in bl))),
                   polish_bad_sign=int(sum(r["polish"]["bad_sign"] for r in rs)),
                   active_share=float(np.median([r["polish"]["n_active"] / (len(r["t"]) * r["n_live"]) for r in rs])),
                   depth=depth_tables(rs))
        row["err"] = {k: row["alpha"][k] - K for k in keys}
        out["profiles"][p] = row
        dF += [abs(d["F"] - d["R"]) for d in ds]
        inter += [d["F_q25"] < d["F"] < d["F_q75"] for d in ds]
    P = out["profiles"]
    m1 = {p: abs(P[p]["err"]["R"]) <= 0.05 for p in P}
    out["verdicts"]["M1"] = dict(rows=m1, n_ok=sum(m1.values()), holds=sum(m1.values()) >= 10)
    out["verdicts"]["M2"] = dict(median_abs=float(np.median(dF)), share_le_0_05=float(np.mean(np.asarray(dF) <= 0.05)),
                                 holds=bool(np.median(dF) <= 0.02 and np.mean(np.asarray(dF) <= 0.05) >= 0.9))
    out["verdicts"]["M3"] = dict(max_abs_err={k: float(max(abs(P[p]["err"][k]) for p in P)) for k in ("U", "E", "Y", "H")},
                                 holds_H=bool(all(abs(P[p]["err"]["H"]) <= 0.10 for p in P)))
    out["verdicts"]["M4"] = dict(intervention_order_share=float(np.mean(inter)), holds=bool(np.mean(inter) >= 0.9),
                                 mu_share_of_Jvv={p: 1 - P[p]["blocks"]["share_C"]["tt"] for p in P})
    out["summary_max_abs_err"] = {k: float(max(abs(P[p]["err"][k]) for p in P)) for k in keys}
    out["summary_max_abs_err_r_le_1"] = {k: float(max(abs(P[p]["err"][k]) for p in R_LE_1 if p in P)) for k in keys}
    OUT.mkdir(exist_ok=True); (OUT / "hessian.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    cols = ["U", "E", "Y", "N", "Z", "H", "Hd", "diagH", "diagHd", "F_fix", "F_amp", "F", "R", "med", "F_q25", "F_q75"]
    print("prof   r     K    " + " ".join(f"{c:>7s}" for c in cols))
    for p, row in P.items():
        print(f"{p:5s} {row['r']:4.2f} {row['K_ref']:.3f} " + " ".join(f"{row['err'][c]:+7.3f}" for c in cols)
              + f"   n={row['n']} conv={row['all_converged']}")
    print("max |err|      " + " ".join(f"{out['summary_max_abs_err'][c]:7.3f}" for c in cols))
    for p, row in P.items():
        b = row["blocks"]
        print(p, "C share tt/tm/mm", {k: round(v, 2) for k, v in b["share_C"].items()},
              "alpha", {k: round(v, 3) for k, v in b["alpha_block"].items()}, "closure", round(b["closure_max"], 4))
    print(json.dumps(out["verdicts"], indent=1))


if __name__ == "__main__":
    main()
