"""Pre-registered mechanical judge for T1.1 (PLAN_T11_ZH.md)."""
import json, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "bosr_t11.json"


def main():
    d = json.loads(RES.read_text(encoding="utf-8"))
    verdict = {}

    # 1. M0_repro_fixed
    pa = d["PartA_M0_repro"]
    q_pz = pa["incoh_pz_pow"]["invQ_noint"]
    q_a = pa["incoh_a"]["invQ_noint"]
    q_b = pa["incoh_b"]["invQ_noint"]
    m0_ok = (abs(q_pz - (-7.80)) <= 0.15) and (abs(q_a - (-1.94)) <= 0.10) and (abs(q_b - (-16.70)) <= 0.10)
    verdict["M0_repro_fixed"] = dict(
        pass_rule=bool(m0_ok), c0sq_allwin=pa["c0sq_allwin"],
        pz_pow=q_pz, a=q_a, b=q_b,
    )

    # 2. T11_S0_Null (synthetic Q=inf with velocity error + multiples + ghosts + near noise)
    s0 = d["PartD_Synthetic"]["S0_Qinf"]
    s0_ab150 = s0["H150"]["AutoBlock"]["invQ_fe"]
    s0_ab200 = s0["H200"]["AutoBlock"]["invQ_fe"]
    s0_fix200 = s0["H200"]["Fixed1D"]["invQ_fe"]
    s0_ok = (abs(s0_ab150) <= 2.5) and (abs(s0_ab200) <= 2.5)
    verdict["T11_S0_Null"] = dict(
        pass_rule=bool(s0_ok),
        Incoh_noint=s0["Incoh"]["invQ_noint"], Incoh_fe=s0["Incoh"]["invQ_fe"],
        Fixed1D_H200_fe=s0_fix200,
        AutoBlock_H100_fe=s0["H100"]["AutoBlock"]["invQ_fe"],
        AutoBlock_H150_fe=s0_ab150, AutoBlock_H200_fe=s0_ab200,
    )

    # 3. T11_S1_Recov (synthetic Q_VSP=(60,90,100))
    s1 = d["PartD_Synthetic"]["S1_QVSP"]
    truth = s1["true_invQ_eff"]
    s1_incoh = s1["Incoh"]["invQ_noint"]
    s1_ab100 = s1["H100"]["AutoBlock"]["invQ_fe"]
    s1_ab150 = s1["H150"]["AutoBlock"]["invQ_fe"]
    s1_ab200 = s1["H200"]["AutoBlock"]["invQ_fe"]
    ci150 = s1["H150"]["AutoBlock"]["ci90_fe"]
    drift_s1 = abs(s1_ab200 - s1_ab100)
    s1_ok = (
        (truth - s1_incoh >= 4.0)
        and (abs(s1_ab150 - truth) / max(abs(truth), 1e-9) <= 0.20)
        and (ci150[0] <= truth <= ci150[1])
        and (drift_s1 <= 2.5)
    )
    verdict["T11_S1_Recov"] = dict(
        pass_rule=bool(s1_ok), truth=truth,
        Incoh_noint=s1_incoh, Incoh_fe=s1["Incoh"]["invQ_fe"],
        Fixed1D_H200_fe=s1["H200"]["Fixed1D"]["invQ_fe"],
        AutoBlock_H100_fe=s1_ab100, AutoBlock_H150_fe=s1_ab150, AutoBlock_H200_fe=s1_ab200,
        ci90_H150=ci150, aperture_drift_200_100=drift_s1,
    )

    # 4. T11_Proc_Ctrl (processed Volve seq079_dline5081 kinematics calibration)
    pb = d["PartB_Proc_Ctrl"]
    p_inc = pb["Incoh"]["invQ_noint"]
    p_inc_fe = pb["Incoh"]["invQ_fe"]
    p_ab100 = pb["H100"]["AutoBlock"]["invQ_fe"]
    p_ab150 = pb["H150"]["AutoBlock"]["invQ_fe"]
    p_ab200 = pb["H200"]["AutoBlock"]["invQ_fe"]
    drift_p = abs(p_ab200 - p_ab100)
    proc_ok = (abs(p_ab150 - p_inc) <= 3.5) and (drift_p <= 3.5)
    verdict["T11_Proc_Ctrl"] = dict(
        pass_rule=bool(proc_ok),
        Incoh_noint=p_inc, Incoh_fe=p_inc_fe,
        Fixed1D_H100_fe=pb["H100"]["Fixed1D"]["invQ_fe"],
        Fixed1D_H200_fe=pb["H200"]["Fixed1D"]["invQ_fe"],
        AutoBlock_H100_fe=p_ab100, AutoBlock_H150_fe=p_ab150, AutoBlock_H200_fe=p_ab200,
        aperture_drift_200_100=drift_p,
    )

    # 5. T11_Raw_TrueQ (raw Volve a or pz_up genuine physical Q vs beam-smearing artifact)
    pc = d["PartC_Raw_Volve"]
    raw_details = {}
    any_raw_ok = False
    for comp in ("a", "b", "pz_up"):
        r100 = pc[comp]["H100"]["AutoBlock"]
        r150 = pc[comp]["H150"]["AutoBlock"]
        r200 = pc[comp]["H200"]["AutoBlock"]
        drift = abs(r200["invQ_fe"] - r100["invQ_fe"])
        ok_c = (
            (r100["invQ_fe"] >= 5.0) and (r100["ci90_fe"][0] > 0.0)
            and (r150["invQ_fe"] >= 5.0) and (r150["ci90_fe"][0] > 0.0)
            and (r150["q_lo_12_30_fe"] > 0.0)
            and (drift <= 5.0)
        )
        raw_details[comp] = dict(
            pass_rule=bool(ok_c),
            Fixed1D_H100_fe=pc[comp]["H100"]["Fixed1D"]["invQ_fe"],
            Fixed1D_H200_fe=pc[comp]["H200"]["Fixed1D"]["invQ_fe"],
            AutoBlock_H75_fe=pc[comp]["H75"]["AutoBlock"]["invQ_fe"],
            AutoBlock_H100_fe=r100["invQ_fe"], ci90_H100=r100["ci90_fe"],
            AutoBlock_H150_fe=r150["invQ_fe"], ci90_H150=r150["ci90_fe"],
            AutoBlock_H200_fe=r200["invQ_fe"], ci90_H200=r200["ci90_fe"],
            H150_lo_12_30_fe=r150["q_lo_12_30_fe"], H150_hi_30_55_fe=r150["q_hi_30_55_fe"],
            H150_c2_over_dS=r150["c2_over_dS_med"],
            H150_alpha_med=r150["alpha_med"],
            aperture_drift_200_100=drift,
        )
        if comp in ("a", "pz_up") and ok_c:
            any_raw_ok = True
    verdict["T11_Raw_TrueQ"] = dict(pass_rule=bool(any_raw_ok), details=raw_details)

    out_path = ROOT / "results" / "bosr_t11_verdict.json"
    out_path.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
