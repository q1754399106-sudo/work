"""Pre-registered mechanical judge for T3 (anchor_norm_phasediag_20261123/PLAN_ZH.md)."""
import json
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
    closures = [r["rel_closure"] for r in recs if np.isfinite(r.get("rel_closure", np.nan))]
    n0_ok = bool(closures and max(closures) <= 1e-3)

    errs_all = [abs(r["alpha_F"] - r["kappa_hat"]) for r in valid_recs if np.isfinite(r["alpha_F"])]
    errs_mu = [abs(r["alpha_F"] - r["kappa_hat"]) for r in valid_recs if r["rule"] in ("R_med", "R_unit") and np.isfinite(r["alpha_F"])]
    by_rule = {}
    for rule in ("R_fix", "R_q25", "R_med", "R_q75", "R_unit"):
        sub = [r for r in recs if r["rule"] == rule]
        vsub = [r for r in valid_recs if r["rule"] == rule]
        by_rule[rule] = dict(
            n_total=len(sub), n_valid=len(vsub),
            n_edge_low=sum(1 for r in sub if r["scans"]["a1.0"]["v_sel"] <= 0.405),
            v_sel_10_med=float(np.median([r["v_sel_10"] for r in sub])) if sub else None,
            kappa_hat_med=float(np.median([r["kappa_hat"] for r in vsub])) if vsub else None,
            alpha_F_med=float(np.median([r["alpha_F"] for r in vsub])) if vsub else None,
            abs_err_med=float(np.median([abs(r["alpha_F"] - r["kappa_hat"]) for r in vsub])) if vsub else None,
            abs_err_max=float(np.max([abs(r["alpha_F"] - r["kappa_hat"]) for r in vsub])) if vsub else None,
        )
    n1_ok = bool(errs_all and float(np.median(errs_all)) <= 0.05 and errs_mu and float(np.median(errs_mu)) <= 0.04)

    # N2_PhaseBoundary
    idx = {(r["profile"], r["rule"], round(r["tw"], 1), r["seed"]): r for r in recs}
    rfix_all = [r for r in recs if r["rule"] == "R_fix"]
    rfix_collapse_frac = float(sum(1 for r in rfix_all if r["scans"]["a1.0"]["v_sel"] <= 0.405) / max(len(rfix_all), 1))

    q_order_ok = True
    q_table = {}
    for prof in ("S", "D", "N3"):
        v25 = idx.get((prof, "R_q25", 1.4, 9931), {}).get("v_sel_10")
        v50 = idx.get((prof, "R_med", 1.4, 9931), {}).get("v_sel_10")
        v75 = idx.get((prof, "R_q75", 1.4, 9931), {}).get("v_sel_10")
        vu  = idx.get((prof, "R_unit", 1.4, 9931), {}).get("v_sel_10")
        q_table[prof] = dict(R_q25=v25, R_med=v50, R_q75=v75, R_unit=vu)
        if None in (v25, v50, v75) or not (v25 < v50 < v75):
            q_order_ok = False

    dv_unit_by_prof = {}
    for prof in ("S", "D", "N3"):
        dvs = []
        for (p, rule, tw, sd), r in idx.items():
            if p == prof and rule == "R_med":
                ru = idx.get((p, "R_unit", tw, sd))
                if ru:
                    dvs.append(ru["v_sel_10"] - r["v_sel_10"])
        dv_unit_by_prof[prof] = dict(med=float(np.median(dvs)) if dvs else None, vals=dvs)

    # N3_DepthScaling (seed=9931, tw = 1.4 -> 2.0 -> 2.5 on R_med & R_unit)
    depth_table = {}
    mono_ok = True
    for rule in ("R_med", "R_unit", "R_fix"):
        depth_table[rule] = {}
        for prof in ("S", "D", "N3"):
            vals = {}
            for tw in (1.4, 2.0, 2.5):
                r = idx.get((prof, rule, tw, 9931))
                vals[f"tw{tw:.1f}"] = None if r is None else dict(
                    v10=r["v_sel_10"], khat=r["kappa_hat"], aF=r["alpha_F"], valid=r["valid"], a_tw=r["a_tw"]
                )
            depth_table[rule][prof] = vals
            if rule in ("R_med", "R_unit"):
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
            med_abs_err_med_unit=float(np.median(errs_mu)) if errs_mu else None,
            by_rule=by_rule,
        ),
        N2_PhaseBoundary=dict(
            rfix_collapse_frac=rfix_collapse_frac,
            quantile_order_pass=bool(q_order_ok),
            quantile_table_tw14_s9931=q_table,
            dv_Runit_minus_Rmed_by_prof=dv_unit_by_prof,
        ),
        N3_DepthScaling=dict(pass_rule=bool(mono_ok), table=depth_table),
    )
    out_path = RES_DIR / "phasediag_verdict.json"
    out_path.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
