"""Pre-registered verdict script for T1.0 (PLAN_ZH.md section 2). Committed before running bosr_diag.py."""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
inp = ROOT / "results" / "bosr_diag.json"
d = json.loads(inp.read_text(encoding="utf-8"))

incoh = d["D0_D2_incoh"]
bosr = d["D3_bosr"]
well = d["D4_well_forward"]
curv = d["D5_curvature"]

# M0: Exact reproduction of field79raw
m0 = (
    abs(incoh["pz_pow"]["invQ_noint"] - (-7.80)) <= 0.10
    and abs(incoh["a"]["invQ_noint"] - (-1.94)) <= 0.10
    and abs(incoh["b"]["invQ_noint"] - (-16.70)) <= 0.10
)

# M-Near: FE intercept lifts 1.2-1.6 s raw PZ_pow above +5.0
m_near = incoh["pz_pow"]["q_12_16_fe"] > 5.0

# M-Mult: on 2.0-2.2 s, multiple beam is bluer than primary (diff <= -0.005 /Hz) AND B-OSR lifts 2.0-2.2 s by >= +8.0 (or > 0)
m_mult_details = {}
m_mult_any = False
for comp in ("a", "pz_up"):
    ss = bosr[comp].get("slowscan_2.0_2.2", {})
    sl_diff = ss.get("diff_mult_minus_prim_med", 0.0)
    q_old = incoh[comp]["q_20_22_noint"]
    q_new = max(bosr[comp]["q_20_22_noint"], bosr[comp]["q_20_22_fe"])
    lift = q_new - q_old
    ok = (sl_diff <= -0.005) and (lift >= 8.0 or q_new > 0.0)
    m_mult_details[comp] = dict(sl_diff_mult_minus_prim=sl_diff, q_old=q_old, q_new=q_new, lift=lift, pass_rule=ok)
    m_mult_any = m_mult_any or ok

# M-Tune: B-OSR still < -5 on 2.0-2.2 s AND well Q=inf gives <= -5 on 2.0-2.2 s
well_qinf_20_22 = {w: well[w]["Q_inf"]["q_20_22"] for w in well}
bosr_20_22_best = max(bosr["a"]["q_20_22_noint"], bosr["a"]["q_20_22_fe"],
                      bosr["pz_up"]["q_20_22_noint"], bosr["pz_up"]["q_20_22_fe"])
m_tune = (bosr_20_22_best < -5.0) and any(v <= -5.0 for v in well_qinf_20_22.values())

# M-Recover: overall B-OSR (noint or fe) on A or PZ_up falls inside VSP [5.6, 17.4] with 90% CI width < 20
m_recover_details = {}
m_recover_any = False
for comp in ("a", "pz_up"):
    for mode, qk, ck in [("noint", "invQ_noint", "ci90_noint"), ("fe", "invQ_fe", "ci90_fe")]:
        qv = bosr[comp][qk]
        ci = bosr[comp][ck]
        width = ci[1] - ci[0]
        ok = (5.6 <= qv <= 17.4) and (width < 20.0)
        m_recover_details[f"{comp}_{mode}"] = dict(q=qv, ci90=ci, width=width, pass_rule=ok)
        m_recover_any = m_recover_any or ok

# Curv-Proc: raw near-offset c2 < +4.0e-4 (30% of +1.322e-3)
c2_raw_up = curv.get("raw_PZ_up", {}).get("c2", float("nan"))
c2_raw_a = curv.get("raw_A_TRID11", {}).get("c2", float("nan"))
curv_proc = (c2_raw_up < 4.0e-4) or (c2_raw_a < 4.0e-4)

verdict = dict(
    M0_repro=m0,
    M_Near=dict(pass_rule=m_near, pz_pow_12_16_noint=incoh["pz_pow"]["q_12_16_noint"],
                pz_pow_12_16_fe=incoh["pz_pow"]["q_12_16_fe"]),
    M_Mult=dict(pass_rule=m_mult_any, details=m_mult_details),
    M_Tune=dict(pass_rule=m_tune, bosr_20_22_best=bosr_20_22_best, well_Qinf_20_22=well_qinf_20_22),
    M_Recover=dict(pass_rule=m_recover_any, details=m_recover_details),
    Curv_Proc=dict(pass_rule=curv_proc, c2_raw_PZ_up=c2_raw_up, c2_raw_A=c2_raw_a,
                   c2_proc_near=curv.get("proc_seq079_near", {}).get("c2"),
                   c2_poststack_ref=curv.get("poststack_archived_ref_c2")),
)
out_path = ROOT / "results" / "bosr_verdict.json"
out_path.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(verdict, indent=2, ensure_ascii=False))
