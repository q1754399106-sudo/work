"""Scores of PLAN_VARPRO_ZH.md -> analysis/varpro_curv.json and a printed table.
VP: alpha = a0 sum num / sum den (stored per dataset as alpha_vp). NF: phi = 0.1 s moving average of n(t) divided by the
atom spectral centroid fc(t). G (descriptive): phi = phiG. E, N reproduced against results/ (archive)."""
import json
from pathlib import Path
import numpy as np
from analyze import REF, R_LE_1, smooth, alpha, uniform   # same windows and reference values as PLAN_ZH.md

HERE = Path(__file__).resolve().parents[1]; RES = HERE / "results_vp"; ARC = HERE / "results"; OUT = HERE / "analysis"


def a_phi(rec, phi):
    a = np.asarray(rec["a"], float); return float(rec["a0"] * np.sum(phi * a) / np.sum(phi * a ** 2))


def main():
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(RES.glob("*.json"))]
    rep = []; groups = {}
    for r in recs:
        arc = json.loads((ARC / f"{r['profile']}_{r['id']}.json").read_text(encoding="utf-8"))
        rep.append(max(abs(alpha(r, "e") - alpha(arc, "e")), abs(alpha(r, "n") - alpha(arc, "n"))))
        nf = a_phi(r, smooth(r["n"]) / np.asarray(r["fc"], float)); g = a_phi(r, np.asarray(r["phiG"], float))
        groups.setdefault(r["profile"], []).append(dict(VP=r["alpha_vp"], NF=nf, G=g, E=alpha(r, "e"), N=alpha(r, "n"), U=uniform(r),
                                                      ratio=r["ratio_1phi1_tr"], conv=r["converged"]))
    out = dict(n=len(recs), reproduction_max_abs=float(max(rep)), reproduction_ok=bool(max(rep) <= 0.002), profiles={}, verdicts={})
    print("datasets %d, reproduction max |d alpha| %.4f" % (len(recs), max(rep)))
    print("%-5s %5s %6s %7s | %7s %7s %7s %7s %7s | ratio" % ("prof", "r", "K", "e_U", "e_VP", "e_NF", "e_G", "e_E", "e_N"))
    for p in REF:
        if p not in groups:
            continue
        r_, aU, K = REF[p]; g = groups[p]; med = {k: float(np.median([d[k] for d in g])) for k in ("VP", "NF", "G", "E", "N", "U")}
        rng = {k: [float(min(d[k] for d in g)), float(max(d[k] for d in g))] for k in ("VP", "NF")}
        e = {k: med[k] - K for k in med}; eU = aU - K
        out["profiles"][p] = dict(r=r_, K=K, alpha_U=aU, e_U=eU, alpha=med, e=e, range=rng, n=len(g),
                                  ratio_1phi1_tr_median=float(np.median([d["ratio"] for d in g])),
                                  not_converged=int(sum(not d["conv"] for d in g)), nf_minus_U=med["NF"] - med["U"])
        print("%-5s %5.2f %6.3f %+7.3f | %+7.3f %+7.3f %+7.3f %+7.3f %+7.3f | %.2f" % (p, r_, K, eU, e["VP"], e["NF"], e["G"], e["E"], e["N"],
                                                                                    out["profiles"][p]["ratio_1phi1_tr_median"]))
    for k in ("VP", "NF"):
        P = out["profiles"]
        a_ok = all(abs(P[p]["e"][k]) <= abs(P[p]["e_U"]) + 0.02 for p in R_LE_1)
        b_ok = abs(P["D"]["e"][k]) <= 0.152
        out["verdicts"][k] = dict(a=a_ok, b=b_ok, a_fail=[p for p in R_LE_1 if abs(P[p]["e"][k]) > abs(P[p]["e_U"]) + 0.02])
        print(k, "a", a_ok, "b", b_ok, "a fails:", out["verdicts"][k]["a_fail"])
    P = out["profiles"]
    out["nf_reduces_to_U_on_r_le_1"] = all(abs(P[p]["nf_minus_U"]) <= 0.02 for p in R_LE_1)
    OUT.mkdir(exist_ok=True); (OUT / "varpro_curv.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
