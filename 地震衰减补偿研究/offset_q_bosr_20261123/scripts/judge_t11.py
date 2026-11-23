"""Pre-registered mechanical judge for T1.1 (PLAN_T11_ZH.md with CW-BOSR)."""
import json
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
    verdict["M0_repro_fixed"] = dict(pass_rule=bool(m0_ok), c0sq_allwin=pa["c0sq_allwin"], pz_pow=q_pz, a=q_a, b=q_b)

    # 2. T11_S0_Null (synthetic Q=inf Earth with velocity error + multiples + ghosts + near noise)
    s0 = d["PartD_Synthetic"]["S0_Qinf"]
    s0_cw36_ni = s0["CW3600"]["invQ_noint"]; s0_cw36_fe = s0["CW3600"]["invQ_fe"]
    s0_cw48_ni = s0["CW4800"]["invQ_noint"]; s0_cw48_fe = s0["CW4800"]["invQ_fe"]
    s0_ok = all(abs(v) <= 3.0 for v in (s0_cw36_ni, s0_cw36_fe, s0_cw48_ni, s0_cw48_fe))
    verdict["T11_S0_Null"] = dict(
        pass_rule=bool(s0_ok),
        Incoh_noint=s0["Incoh"]["invQ_noint"], Incoh_fe=s0["Incoh"]["invQ_fe"],
        Fix_H200_fe=s0["H200"]["Fixed1D"]["invQ_fe"],
        CW3600_noint=s0_cw36_ni, CW3600_fe=s0_cw36_fe,
        CW4800_noint=s0_cw48_ni, CW4800_fe=s0_cw48_fe,
    )

    # 3. T11_S1_Recov (synthetic Q_VSP=(60,90,100))
    s1 = d["PartD_Synthetic"]["S1_QVSP"]
    truth = s1["true_invQ_eff"]
    s1_inc = s1["Incoh"]["invQ_noint"]
    s1_cw36 = s1["CW3600"]; s1_cw48 = s1["CW4800"]
    ci36 = s1_cw36["ci90_noint"]
    drift_s1 = abs(s1_cw48["invQ_noint"] - s1_cw36["invQ_noint"])
    s1_ok = (
        (truth - s1_inc >= 10.0)
        and (abs(s1_cw36["invQ_noint"] - truth) / max(abs(truth), 1e-9) <= 0.20)
        and (ci36[0] <= truth <= ci36[1])
        and (drift_s1 <= 3.5)
    )
    verdict["T11_S1_Recov"] = dict(
        pass_rule=bool(s1_ok), truth=truth,
        Incoh_noint=s1_inc, Incoh_fe=s1["Incoh"]["invQ_fe"],
        Fix_H200_noint=s1["H200"]["Fixed1D"]["invQ_noint"], Fix_H200_fe=s1["H200"]["Fixed1D"]["invQ_fe"],
        CW3600_noint=s1_cw36["invQ_noint"], CW3600_ci90=ci36, CW3600_fe=s1_cw36["invQ_fe"],
        CW4800_noint=s1_cw48["invQ_noint"], CW4800_fe=s1_cw48["invQ_fe"],
        scale_drift_4800_3600=drift_s1,
    )

    # 4. T11_Proc_Ctrl (processed Volve seq079_dline5081 scale invariance)
    pb = d["PartB_Proc_Ctrl"]
    p_drift = abs(pb["CW4800"]["invQ_fe"] - pb["CW3600"]["invQ_fe"])
    p_amed = pb["H200"]["AutoBlock"]["alpha_med"]
    proc_ok = (p_drift <= 4.0) and (0.92 <= p_amed <= 1.08)
    verdict["T11_Proc_Ctrl"] = dict(
        pass_rule=bool(proc_ok),
        Incoh_noint=pb["Incoh"]["invQ_noint"], Incoh_fe=pb["Incoh"]["invQ_fe"],
        Fix_H100_fe=pb["H100"]["Fixed1D"]["invQ_fe"], Fix_H200_fe=pb["H200"]["Fixed1D"]["invQ_fe"],
        CW3600_noint=pb["CW3600"]["invQ_noint"], CW3600_fe=pb["CW3600"]["invQ_fe"],
        CW4800_noint=pb["CW4800"]["invQ_noint"], CW4800_fe=pb["CW4800"]["invQ_fe"],
        H200_Auto_alpha_med=p_amed, scale_drift_4800_3600=p_drift,
    )

    # 5. T11_Raw_TrueQ (raw Volve CW-BOSR artifact-free check)
    pc = d["PartC_Raw_Volve"]
    raw_details = {}
    any_ok = False
    for comp in ("a", "b", "pz_up"):
        c36 = pc[comp]["CW3600"]
        c48 = pc[comp]["CW4800"]
        drift = abs(c48["invQ_fe"] - c36["invQ_fe"])
        ok_c = (
            (c36["invQ_fe"] >= 5.0) and (c36["ci90_fe"][0] > 0.0)
            and (c48["invQ_fe"] >= 5.0) and (c48["ci90_fe"][0] > 0.0)
            and (c36["q_lo_12_30_fe"] > 0.0)
            and (drift <= 4.0)
        )
        raw_details[comp] = dict(
            pass_rule=bool(ok_c),
            Fix_H100_fe=pc[comp]["H100"]["Fixed1D"]["invQ_fe"],
            Fix_H150_fe=pc[comp]["H150"]["Fixed1D"]["invQ_fe"],
            Fix_H200_fe=pc[comp]["H200"]["Fixed1D"]["invQ_fe"],
            CW3600_noint=c36["invQ_noint"], CW3600_ci90_noint=c36["ci90_noint"],
            CW3600_fe=c36["invQ_fe"], CW3600_ci90_fe=c36["ci90_fe"],
            CW3600_lo_hi_fe=[c36["q_lo_12_30_fe"], c36["q_hi_30_55_fe"]],
            CW3600_1216_2022_fe=[c36["q_12_16_fe"], c36["q_20_22_fe"]],
            CW4800_noint=c48["invQ_noint"], CW4800_fe=c48["invQ_fe"], CW4800_ci90_fe=c48["ci90_fe"],
            CW4800_lo_hi_fe=[c48["q_lo_12_30_fe"], c48["q_hi_30_55_fe"]],
            scale_drift_4800_3600=drift,
        )
        if comp in ("a", "pz_up") and ok_c:
            any_ok = True
    verdict["T11_Raw_TrueQ"] = dict(pass_rule=bool(any_ok), details=raw_details)

    out_path = ROOT / "results" / "bosr_t11_verdict.json"
    out_path.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
