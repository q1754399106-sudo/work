"""Pre-registered automated mechanical judge for uq_rtm_unified_20261124 (PLAN_ZH.md)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = RES / "unified_verdict.json"

S4_F25 = Path(r"D:\ChatGPT\研究\地震衰减补偿研究\uq_rtm_20261117\analysis\stage4_field_fref25.json")
S4H_JSON = Path(r"D:\ChatGPT\研究\.claude\worktrees\study-research-docs-plan-c4d4c6\地震衰减补偿研究\uq_rtm_20261117\analysis\stage4h_field.json")


def _arr(x):
    return [float(v) for v in x] if isinstance(x, list) else [float(x)]


def _load(p):
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def judge_tier1():
    u_models = {m: _load(RES / f"synth_{m}.json") for m in ["U1", "U2", "U3", "U4", "U5"]}
    out = {}

    # G1: Theorem 3 B=300 in DSA domain across U1-U3 (l30)
    if all(u_models[m] is not None for m in ["U1", "U2", "U3"]):
        raw_b1 = {m: float(u_models[m]["noise_acc"]["RAW_B1_log2err"]) for m in ["U1", "U2", "U3"]}
        dsa_b300 = {m: float(u_models[m]["noise_acc"]["DSA_B300_log2err"]) for m in ["U1", "U2", "U3"]}
        raw_b300 = {m: float(u_models[m]["noise_acc"]["RAW_B300_log2err"]) for m in ["U1", "U2", "U3"]}
        dsa_b1 = {m: float(u_models[m]["noise_acc"]["DSA_B1_log2err"]) for m in ["U1", "U2", "U3"]}
        g1_pass = all(v <= -2.0 for v in raw_b1.values()) and all(abs(v) <= 1.0 for v in dsa_b300.values())
        out["G1_Theorem3_B300_In_DSA"] = {
            "PASS": g1_pass,
            "RAW_B1_log2err": raw_b1,
            "RAW_B300_log2err": raw_b300,
            "DSA_B1_log2err": dsa_b1,
            "DSA_B300_log2err": dsa_b300,
        }

        # G2: Confluence beats every partial component on U1-U3 (W_deep)
        d_vs_raw_b1, d_vs_raw_b300, d_vs_dsa_b1, d_vs_f100, d_rob_m05 = {}, {}, {}, {}, {}
        for m in ["U1", "U2", "U3"]:
            wd = u_models[m]["W_deep"]
            hero = float(wd["DSA_UQW1b_25"]["E1_mean"])
            d_vs_raw_b1[m] = hero - float(wd["RAW_UQW1_25"]["E1_mean"])
            d_vs_raw_b300[m] = hero - float(wd["RAW_UQW1b_25"]["E1_mean"])
            d_vs_dsa_b1[m] = hero - float(wd["DSA_UQW1_25"]["E1_mean"])
            d_vs_f100[m] = hero - float(wd["RAW_UQW1b_100"]["E1_mean"])
            h05 = _arr(wd["DSA_UQW1b_25"]["E1"])[0]
            c05 = max(_arr(wd["DSA_IQ70_25"]["E1"])[0], _arr(wd["DSA_W1b_25"]["E1"])[0])
            d_rob_m05[m] = h05 - c05

        c1 = (sum(d_vs_raw_b1.values()) / 3.0 >= 0.025) and all(v > 0 for v in d_vs_raw_b1.values())
        c2 = (sum(d_vs_raw_b300.values()) / 3.0 >= 0.015) and all(v > 0 for v in d_vs_raw_b300.values())
        c3 = (sum(d_vs_dsa_b1.values()) / 3.0 >= 0.010) and all(v > 0 for v in d_vs_dsa_b1.values())
        c4 = (sum(d_vs_f100.values()) / 3.0 >= 0.150) and all(v > 0 for v in d_vs_f100.values())
        c5 = (sum(d_rob_m05.values()) / 3.0 >= 0.030) and all(v > 0 for v in d_rob_m05.values())
        out["G2_Synth_Confluence_Beats_Components"] = {
            "PASS": c1 and c2 and c3 and c4 and c5,
            "sub_checks": {"vs_RAW_B1": c1, "vs_RAW_B300": c2, "vs_DSA_B1": c3, "vs_fref100": c4, "m05_robust": c5},
            "delta_vs_RAW_UQW1_25": d_vs_raw_b1,
            "delta_vs_RAW_UQW1b_25": d_vs_raw_b300,
            "delta_vs_DSA_UQW1_25": d_vs_dsa_b1,
            "delta_vs_RAW_UQW1b_100": d_vs_f100,
            "delta_m05_vs_max_IQ70_W1b": d_rob_m05,
        }

    # G3: U4 Volve full-pipeline twin
    if u_models["U4"] is not None:
        M = u_models["U4"]["methods"]
        r0 = float(M["RAW_R0"]["mean"])
        raw_uq = _arr(M["RAW_UQW1_25"]["tie"])
        hero = _arr(M["DSA_UQW1b_PAnc_25"]["tie"])
        w1b = _arr(M["DSA_W1b_PAnc_25"]["tie"])
        iq70 = _arr(M["DSA_IQ70_PAnc_25"]["tie"])
        repro_drop = (sum(raw_uq) / len(raw_uq)) < r0
        beats_raw = ((sum(hero) / len(hero) - sum(raw_uq) / len(raw_uq)) >= 0.050) and ((hero[2] - raw_uq[2]) >= 0.050)
        beats_r0 = (sum(hero) / len(hero) >= r0) and (hero[2] >= r0)
        beats_m05 = (hero[0] > w1b[0]) and (hero[0] > iq70[0])
        out["G3_Synth_VolveTwin_U4_FullConfluence"] = {
            "PASS": repro_drop and beats_raw and beats_r0 and beats_m05,
            "sub_checks": {"repro_drop": repro_drop, "beats_raw_by_005": beats_raw, "beats_R0": beats_r0, "m05_beats_det": beats_m05},
            "RAW_R0": r0,
            "RAW_UQW1_25_mean": sum(raw_uq) / len(raw_uq),
            "DSA_UQW1b_PAnc_25_mean": sum(hero) / len(hero),
            "DSA_UQW1b_PAnc_25_m05_m10_m20": [hero[0], hero[2], hero[4]],
            "DSA_W1b_PAnc_25_m05": w1b[0],
            "DSA_IQ70_PAnc_25_m05": iq70[0],
        }

    # G7: U5 High-SNR + velocity error safety
    if u_models["U5"] is not None:
        wd = u_models["U5"]["W_deep"]
        de1 = float(wd["DSA_UQW1b_25"]["E1_mean"]) - float(wd["RAW_UQW1_25"]["E1_mean"])
        dsnr = float(wd["DSA_UQW1b_25"]["SNR_img_dB_mean"]) - float(wd["RAW_UQW1_25"]["SNR_img_dB_mean"])
        out["G7_VelocityError_And_HighSNR_Safety"] = {
            "PASS": (de1 >= -0.015) and (dsnr >= 0.5),
            "delta_E1_mean": de1,
            "delta_SNR_dB_mean": dsnr,
        }
    return out


def judge_tier2():
    d = _load(RES / "volve25_unified.json")
    s4 = _load(S4_F25)
    if d is None or s4 is None:
        return {}
    M = d["methods"]
    sl = d["res_slope"]
    s4_idx = [0, 2, 4, 6, 8]
    s4_r0 = float(s4["R0"]["tie_F5"])
    s4_uq25 = [float(s4["UQW1"]["tie_F5"][i]) for i in s4_idx]
    s4_uq25_mean = sum(s4_uq25) / len(s4_uq25)

    r0_new = _arr(M["RAW_R0"]["F5_R0band"])[0]
    r0_br = _arr(M["RAW_R0"]["F5_Broad"])[0]
    hero_r0 = _arr(M["DSA_UQW1b_PAnc_25"]["F5_R0band"])
    hero_br = _arr(M["DSA_UQW1b_PAnc_25"]["F5_Broad"])
    amp_r0 = _arr(M["DSA_UQW1b_PAnc_Amp"]["F5_R0band"])
    amp_br = _arr(M["DSA_UQW1b_PAnc_Amp"]["F5_Broad"])

    sm_raw = _arr(M["RAW_R0"]["Smile_Frac"])[0]
    sm_hero = _arr(M["DSA_UQW1b_PAnc_25"]["Smile_Frac"])[2]
    sl_hero = abs(float(sl["DSA_UQW1b_PAnc_25"]))
    sl_eq = abs(float(sl["DSA_R0_Eq"]))
    sl_un = abs(float(sl["RAW_UQW1b_Unanc_25"]))

    c_repro = abs(r0_new - s4_r0) < 1e-4
    c_smile = sm_hero <= 0.25 * sm_raw
    c_slope = (sl_hero <= 0.5 * sl_eq) and (sl_hero <= 0.5 * sl_un)
    c_beat_s4 = (sum(hero_r0) / len(hero_r0) - s4_uq25_mean) >= 0.100
    c_beat_br = (sum(hero_br) / len(hero_br) > r0_br) or (hero_br[2] > r0_br)
    return {
        "G4_Volve25_Confluence_Closure": {
            "PASS": c_repro and c_smile and c_slope and c_beat_s4 and c_beat_br,
            "sub_checks": {
                "repro_R0": c_repro,
                "smile_cut_75pct": c_smile,
                "slope_cut_50pct": c_slope,
                "beats_S4_UQW1_25_by_010": c_beat_s4,
                "beats_RAW_R0_Broad": c_beat_br,
            },
            "RAW_R0_diff": abs(r0_new - s4_r0),
            "Smile_Frac_RAW_to_Hero": [sm_raw, sm_hero],
            "abs_slope_Hero_Eq_Unanc": [sl_hero, sl_eq, sl_un],
            "S4_UQW1_25_mean": s4_uq25_mean,
            "Hero25_R0band_mean": sum(hero_r0) / len(hero_r0),
            "Hero25_Broad_mean": sum(hero_br) / len(hero_br),
            "Amp_R0band_mean": sum(amp_r0) / len(amp_r0),
            "Amp_Broad_mean": sum(amp_br) / len(amp_br),
            "RAW_R0_Broad": r0_br,
        }
    }


def judge_tier3():
    d = _load(RES / "volve45_unified.json")
    s4h = _load(S4H_JSON)
    if d is None or s4h is None:
        return {}
    M = d["methods"]
    sl = d["res_slope"]

    s4h_r0 = float(s4h["R0"]["tie_fix_F5"])
    s4h_uqb_m1 = float(s4h["UQW1b"]["tie_fix_F5"][2])

    r0_new = _arr(M["RAW_R0"]["F5_fix"])[0]
    uqb_new = _arr(M["RAW_UQW1b_Unanc_100"]["F5_fix"])[2]
    hero25 = _arr(M["DSA_UQW1b_PAnc_25"]["F5_fix"])
    w1b25 = _arr(M["DSA_W1b_PAnc_25"]["F5_fix"])
    iq7025 = _arr(M["DSA_IQ70_PAnc_25"]["F5_fix"])
    amp = _arr(M["DSA_UQW1b_PAnc_Amp"]["F5_fix"])
    h100 = _arr(M["DSA_UQW1b_PAnc_100"]["F5_fix"])

    sm_s4h = _arr(M["RAW_UQW1b_Unanc_100"]["Smile_Frac"])[2]
    sm_hero = _arr(M["DSA_UQW1b_PAnc_25"]["Smile_Frac"])[2]
    sl_hero = abs(float(sl["DSA_UQW1b_PAnc_25"]))
    sl_eq = abs(float(sl["DSA_R0_Eq"]))
    sl_un = abs(float(sl["RAW_UQW1b_Unanc_100"]))

    c_repro = (abs(r0_new - s4h_r0) < 1e-3) and (abs(uqb_new - s4h_uqb_m1) < 1e-3)
    c_smile = sm_hero <= 0.30 * sm_s4h
    c_slope = (sl_hero <= 0.5 * sl_eq) and (sl_hero <= 0.5 * sl_un)
    g5 = {
        "PASS": c_repro and c_smile and c_slope,
        "sub_checks": {"repro_S4H": c_repro, "smile_cut_70pct": c_smile, "slope_cut_50pct": c_slope},
        "repro_R0_abs_diff": abs(r0_new - s4h_r0),
        "repro_UQW1b_m1_abs_diff": abs(uqb_new - s4h_uqb_m1),
        "Smile_Frac_S4H_UQW1b_to_Hero25": [sm_s4h, sm_hero],
        "abs_slope_Hero25_Eq_Unanc100": [sl_hero, sl_eq, sl_un],
    }

    c_all_m = all(v >= r0_new for v in hero25)
    c_m1_03 = (hero25[2] - r0_new) >= 0.030
    c_m06_det = (hero25[0] > w1b25[0]) and (hero25[0] > iq7025[0])
    c_no_col = hero25[0] >= hero25[4] - 0.030
    g6 = {
        "PASS": c_all_m and c_m1_03 and c_m06_det and c_no_col,
        "sub_checks": {
            "all_5m_ge_R0": c_all_m,
            "m1_beats_R0_by_003": c_m1_03,
            "m06_beats_W1b_IQ70": c_m06_det,
            "no_overcomp_collapse": c_no_col,
        },
        "RAW_R0_fix": r0_new,
        "Stage4H_RAW_UQW1b_100_m1": uqb_new,
        "Hero25_F5_fix": hero25,
        "Hero25_F5_fix_mean": sum(hero25) / len(hero25),
        "W1b25_m06": w1b25[0],
        "IQ70_25_m06": iq7025[0],
        "Companion_Amp_F5_fix": amp,
        "Companion_100_F5_fix": h100,
    }
    return {
        "G5_Volve45_Smile_And_Spectral_Closure": g5,
        "G6_Volve45_Broadband_WellTie_And_Robustness": g6,
    }


def main():
    gates = {}
    gates.update(judge_tier1())
    gates.update(judge_tier2())
    gates.update(judge_tier3())
    n_pass = sum(1 for v in gates.values() if v.get("PASS"))
    verdict = {
        "gates": gates,
        "overall_pass_count": n_pass,
        "overall_total": len(gates),
    }
    OUT.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
