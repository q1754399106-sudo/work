"""Pre-registered mechanical judge for T3 (anchor_norm_phasediag_20261123/PLAN_ZH.md)."""
import json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES_DIR = ROOT / "results"


def main():
    files = sorted(RES_DIR.glob("*.json"))
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in files if p.name != "phasediag_verdict.json"]
    if not recs:
        raise SystemExit("No result JSON files found.")

    valid_recs = [r for r in recs if r.get("valid")]
    # N0_Closure
    closures = [r["rel_closure"] for r in valid_recs]
    n0_ok = bool(closures and max(closures) <= 1e-3)

    # N1_ClosedForm
    errs_all = [abs(r["alpha_F"] - r["kappa_hat"]) for r in valid_recs if np.isfinite(r["alpha_F"])]
    errs_mf = [abs(r["alpha_F"] - r["kappa_hat"]) for r in valid_recs if r["rule"] in ("R_med", "R_fix") and np.isfinite(r["alpha_F"])]
    by_rule = {}
    for rule in ("R_med", "R_fix", "R_unit"):
        sub = [r for r in recs if r["rule"] == rule]
        vsub = [r for r in valid_recs if r["rule"] == rule]
        by_rule[rule] = dict(
            n_total=len(sub), n_valid=len(vsub),
            v_sel_10_med=float(np.median([r["v_sel_10"] for r in sub])) if sub else None,
            kappa_hat_med=float(np.median([r["kappa_hat"] for r in vsub])) if vsub else None,
            alpha_F_med=float(np.median([r["alpha_F"] for r in vsub])) if vsub else None,
            abs_err_med=float(np.median([abs(r["alpha_F"] - r["kappa_hat"]) for r in vsub])) if vsub else None,
            abs_err_max=float(np.max([abs(r["alpha_F"] - r["kappa_hat"]) for r in vsub])) if vsub else None,
        )
    n1_ok = bool(errs_all and float(np.median(errs_all)) <= 0.05 and errs_mf and float(np.median(errs_mf)) <= 0.03)

    # N2_RuleShift (paired comparisons at matching (prof, tw, seed))
    idx = {(r["profile"], r["rule"], round(r["tw"], 1), r["seed"]): r for r in recs}
    d_fix_med_k = []
    d_unit_med_v = []
    d_unit_med_k = []
    for (prof, rule, tw, seed), r in idx.items():
        if rule != "R_med":
            continue
        rf = idx.get((prof, "R_fix", tw, seed))
        ru = idx.get((prof, "R_unit", tw, seed))
        if rf and r["valid"] and rf["valid"]:
            d_fix_med_k.append(rf["kappa_hat"] - r["kappa_hat"])
        if ru:
            d_unit_med_v.append(ru["v_sel_10"] - r["v_sel_10"])
            if r["valid"] and ru["valid"]:
                d_unit_med_k.append(ru["kappa_hat"] - r["kappa_hat"])

    # N3_DepthScaling (seed=9931, tw = 1.4 -> 2.0 -> 2.5)
    depth_table = {}
    mono_ok = True
    for rule in ("R_med", "R_fix", "R_unit"):
        depth_table[rule] = {}
        for prof in ("S", "D", "N3"):
            vals = {}
            for tw in (1.4, 2.0, 2.5):
                r = idx.get((prof, rule, tw, 9931))
                vals[f"tw{tw:.1f}"] = None if r is None else dict(
                    v10=r["v_sel_10"], khat=r["kappa_hat"], aF=r["alpha_F"], valid=r["valid"], a_tw=r["a_tw"]
                )
            depth_table[rule][prof] = vals
            if rule in ("R_med", "R_fix"):
                k14 = vals["tw1.4"]["khat"] if vals["tw1.4"] else None
                k20 = vals["tw2.0"]["khat"] if vals["tw2.0"] else None
                k25 = vals["tw2.5"]["khat"] if vals["tw2.5"] else None
                if None in (k14, k20, k25) or not (k25 > k20 > k14):
                    mono_ok = False

    verdict = dict(
        n_total=len(recs), n_valid=len(valid_recs),
        N0_Closure=dict(pass_rule=n0_ok, max_rel_closure=float(max(closures)) if closures else None),
        N1_ClosedForm=dict(
            pass_rule=n1_ok,
            med_abs_err_all=float(np.median(errs_all)) if errs_all else None,
            med_abs_err_med_fix=float(np.median(errs_mf)) if errs_mf else None,
            by_rule=by_rule,
        ),
        N2_RuleShift=dict(
            d_kappa_Rfix_minus_Rmed_med=float(np.median(d_fix_med_k)) if d_fix_med_k else None,
            d_kappa_Rfix_minus_Rmed_vals=d_fix_med_k,
            d_vsel10_Runit_minus_Rmed_med=float(np.median(d_unit_med_v)) if d_unit_med_v else None,
            d_vsel10_Runit_minus_Rmed_vals=d_unit_med_v,
            d_kappa_Runit_minus_Rmed_med=float(np.median(d_unit_med_k)) if d_unit_med_k else None,
        ),
        N3_DepthScaling=dict(pass_rule=bool(mono_ok), table=depth_table),
    )
    out_path = RES_DIR / "phasediag_verdict.json"
    out_path.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
