"""Pre-registered automated judge for uq_rtm_volve_closure_20261123 (PLAN_ZH.md)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES_PATH = ROOT / "results" / "volve_closure_results.json"
STAGE4_AMP_PATH = Path(r"D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\analysis\stage4_field_amp.json")
OUT_PATH = ROOT / "results" / "volve_closure_verdict.json"


def _arr(x):
    return x if isinstance(x, list) else [x]


def main():
    D = json.loads(RES_PATH.read_text(encoding="utf-8"))
    S4 = json.loads(STAGE4_AMP_PATH.read_text(encoding="utf-8"))

    # Stage 4 archived UQW1_amp at m = [0.5, 0.7, 1.0, 1.4, 2.0] (indices 0, 2, 4, 6, 8)
    s4_idx = [0, 2, 4, 6, 8]
    s4_r0 = float(S4["R0"]["tie_F5"])
    s4_uq = [float(S4["UQW1"]["tie_F5"][i]) for i in s4_idx]
    s4_uq_mean = sum(s4_uq) / len(s4_uq)

    M = D["methods"]
    raw_r0_f5 = float(_arr(M["RAW_R0"]["F5_R0band"])[0])
    raw_r0_br = float(_arr(M["RAW_R0"]["F5_Broad"])[0])
    dsa_r0_eq_f5 = float(_arr(M["DSA_R0_Eq"]["F5_R0band"])[0])
    dsa_r0_eq_br = float(_arr(M["DSA_R0_Eq"]["F5_Broad"])[0])

    hero_r0 = [float(v) for v in _arr(M["DSA_UQW1_PAnc"]["F5_R0band"])]
    hero_br = [float(v) for v in _arr(M["DSA_UQW1_PAnc"]["F5_Broad"])]
    hero_r0_mean = sum(hero_r0) / len(hero_r0)
    hero_br_mean = sum(hero_br) / len(hero_br)

    iq_r0 = [float(v) for v in _arr(M["DSA_IQ70_PAnc"]["F5_R0band"])]
    w1_r0 = [float(v) for v in _arr(M["DSA_W1_PAnc"]["F5_R0band"])]

    # Gate V1: Beats Stage 4 archived UQW1_amp by >= +0.050 on mean, m=0.5, m=1.0
    v1_d_mean = hero_r0_mean - s4_uq_mean
    v1_d_m05 = hero_r0[0] - s4_uq[0]
    v1_d_m10 = hero_r0[2] - s4_uq[2]
    v1_pass = (v1_d_mean >= 0.050) and (v1_d_m05 >= 0.050) and (v1_d_m10 >= 0.050)

    # Gate V2: Beats equalized uncompensated baseline DSA_R0_Eq (mean & m=1.0)
    v2_pass = (
        (hero_br_mean > dsa_r0_eq_br)
        and (hero_br[2] > dsa_r0_eq_br)
        and (hero_r0_mean >= dsa_r0_eq_f5)
        and (hero_r0[2] >= dsa_r0_eq_f5)
    )

    # Gate V3: Ultimate closure — beats original Stage-4 RAW_R0 on BOTH R0-band & broadband (at m=1.0 or mean)
    v3_m10 = (hero_r0[2] >= raw_r0_f5) and (hero_br[2] >= raw_r0_br)
    v3_mean = (hero_r0_mean >= raw_r0_f5) and (hero_br_mean >= raw_r0_br)
    v3_pass = v3_m10 or v3_mean

    # Gate V4: Over-compensation robustness at VSP m=0.5
    v4_beats_iq = hero_r0[0] > iq_r0[0]
    v4_beats_w1 = hero_r0[0] > w1_r0[0]
    v4_no_collapse = hero_r0[0] >= hero_r0[4] - 0.030
    v4_pass = v4_beats_iq and v4_beats_w1 and v4_no_collapse

    # Gate V5: Seismic physics closure (spectral flatness + overburden smile suppression)
    sl = D["res_slope"]
    sl_hero = abs(float(sl["DSA_UQW1_PAnc_m1"]))
    sl_eq = abs(float(sl["DSA_R0_Eq"]))
    sl_un = abs(float(sl["RAW_UQW1_Unanc_m1"]))
    v5_spec = (sl_hero <= 0.5 * sl_eq) and (sl_hero <= 0.5 * sl_un)

    sm_raw = float(_arr(M["RAW_R0"]["Smile_Frac"])[0])
    sm_dsa = float(_arr(M["DSA_R0"]["Smile_Frac"])[0])
    v5_smile = sm_dsa < sm_raw
    v5_pass = v5_spec and v5_smile

    verdict = {
        "reproducibility_RAW_R0_vs_Stage4": {
            "RAW_R0_F5_new": raw_r0_f5,
            "RAW_R0_F5_stage4": s4_r0,
            "abs_diff": abs(raw_r0_f5 - s4_r0),
        },
        "gates": {
            "V1_Beats_Stage4_UQW1": {
                "PASS": v1_pass,
                "hero_R0band_mean": hero_r0_mean,
                "stage4_UQW1_amp_mean": s4_uq_mean,
                "delta_mean": v1_d_mean,
                "delta_m05": v1_d_m05,
                "delta_m10": v1_d_m10,
            },
            "V2_Beats_Uncomp_Equalized": {
                "PASS": v2_pass,
                "DSA_R0_Eq_R0band": dsa_r0_eq_f5,
                "hero_R0band_mean": hero_r0_mean,
                "hero_R0band_m10": hero_r0[2],
                "DSA_R0_Eq_Broad": dsa_r0_eq_br,
                "hero_Broad_mean": hero_br_mean,
                "hero_Broad_m10": hero_br[2],
            },
            "V3_Beats_Original_RAW_R0": {
                "PASS": v3_pass,
                "RAW_R0_R0band": raw_r0_f5,
                "hero_R0band_mean": hero_r0_mean,
                "hero_R0band_m10": hero_r0[2],
                "RAW_R0_Broad": raw_r0_br,
                "hero_Broad_mean": hero_br_mean,
                "hero_Broad_m10": hero_br[2],
                "pass_at_m10": v3_m10,
                "pass_at_mean": v3_mean,
            },
            "V4_Overcomp_Robustness_VSP": {
                "PASS": v4_pass,
                "m05_DSA_UQW1_PAnc": hero_r0[0],
                "m05_DSA_W1_PAnc": w1_r0[0],
                "m05_DSA_IQ70_PAnc": iq_r0[0],
                "m20_DSA_UQW1_PAnc": hero_r0[4],
                "m05_minus_m20": hero_r0[0] - hero_r0[4],
            },
            "V5_Seismic_Physics_Closure": {
                "PASS": v5_pass,
                "abs_slope_hero_m1": sl_hero,
                "abs_slope_DSA_R0_Eq": sl_eq,
                "abs_slope_RAW_UQW1_Unanc_m1": sl_un,
                "spec_flat_pass": v5_spec,
                "smile_frac_RAW_R0": sm_raw,
                "smile_frac_DSA_R0": sm_dsa,
                "smile_pass": v5_smile,
            },
        },
        "overall_pass_count": sum([v1_pass, v2_pass, v3_pass, v4_pass, v5_pass]),
        "overall_total": 5,
    }

    OUT_PATH.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
