#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
judge_xmir_qrtm.py
==================
Pre-registered frozen mechanical verdict script for Gates X0-X7
(uq_rtm_xmir_qrtm_20261126). NEVER modify thresholds after pre-registration freeze.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from typing import Any, Dict


def load_j(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    res_dir = os.path.join(root, "results")
    log_dir = os.path.join(root, "logs")
    report_path = os.path.join(root, "RESULTS_ZH.md")

    verdict: Dict[str, Any] = {"gates": {}, "summary": {}}

    # X0: Integrity Gate (Gate R0 / X0 zero-exit audit)
    x0_pass = False
    if os.path.isfile(report_path):
        cmd = [
            sys.executable,
            os.path.join(here, "audit_report_integrity.py"),
            "--report", report_path,
            "--results-dir", res_dir,
            "--scripts-dir", here,
            "--log-out", os.path.join(log_dir, "integrity_audit.log"),
        ]
        rc = subprocess.run(cmd, capture_output=True, text=True).returncode
        x0_pass = (rc == 0)
    verdict["gates"]["X0_Integrity_Gate"] = {"pass": bool(x0_pass), "report_exists": os.path.isfile(report_path)}

    # X1: Analytical Circular-Arc Bent-Ray xmir Proof + Volve Audit (WP1)
    wp1_p = os.path.join(res_dir, "wp1_physics_verify.json")
    if os.path.isfile(wp1_p):
        w1 = load_j(wp1_p)
        t1a = w1["test1a_analytical_ray"]
        reds = {k: float(t1a[k]["error_reduction_pct"]) for k in ("dip_m20", "dip_m10", "dip_p10", "dip_p20")}
        x1_ok = all(v >= 80.0 for v in reds.values()) and ("test1c_volve_xmir_audit" in w1)
        verdict["gates"]["X1_Analytical_BentRay_xmir"] = {
            "pass": bool(x1_ok),
            "dip_error_reduction_pct": reds,
            "volve_res_med_dxmir_m": float(w1["test1c_volve_xmir_audit"]["res_med_dxmir_m"]),
            "volve_res_p90_dxmir_m": float(w1["test1c_volve_xmir_audit"]["res_p90_dxmir_m"]),
            "volve_res_gate_flip_pct": float(w1["test1c_volve_xmir_audit"]["res_gate_flip_pct"]),
        }
    else:
        verdict["gates"]["X1_Analytical_BentRay_xmir"] = {"pass": False, "status": "PENDING_WP1"}

    # Load S1-S4
    s_names = ["S1_Chimney_Base", "S2_StrongGrad_SteepDip", "S3_SharpQ_Multipath", "S4_NonScalar_SpatialQ"]
    s_data: Dict[str, Dict[str, Any]] = {}
    for sn in s_names:
        sp = os.path.join(res_dir, f"synth_{sn}.json")
        if os.path.isfile(sp):
            s_data[sn] = load_j(sp)

    # X2: True Q-RTM No-Strawman Proof on Noise-Free m=1 (S1-S3)
    if all(k in s_data for k in s_names[:3]):
        nf_e1_gains = {}
        nf_eamp_reds = {}
        for sn in s_names[:3]:
            nf = s_data[sn]["NoiseFree_m1"]
            g = float(nf["WE_QRTM_Unreg"]["E1"]) - float(nf["RAW_R0"]["E1"])
            r = 1.0 - float(nf["WE_QRTM_Unreg"]["E_amp"]) / max(float(nf["RAW_R0"]["E_amp"]), 1e-12)
            nf_e1_gains[sn] = g
            nf_eamp_reds[sn] = r
        x2_ok = all(v >= 0.04 for v in nf_e1_gains.values()) and all(v >= 0.60 for v in nf_eamp_reds.values())
        verdict["gates"]["X2_TrueQRTM_NoStrawman_Proof"] = {
            "pass": bool(x2_ok),
            "NF_m1_WE_QRTM_Unreg_E1_gain_vs_R0": nf_e1_gains,
            "NF_m1_WE_QRTM_Unreg_Eamp_red_vs_R0": nf_eamp_reds,
        }
    else:
        verdict["gates"]["X2_TrueQRTM_NoStrawman_Proof"] = {"pass": False, "status": "PENDING_WP2"}

    # X3: Bent-Ray xmir_bent Imaging Benefit on Strong-Gradient Steep-Dip Asymmetric-Cable S2
    if "S2_StrongGrad_SteepDip" in s_data:
        s2 = s_data["S2_StrongGrad_SteepDip"]["W_deep"]
        e1_b = float(s2["DSA_Bent_UQW1b_25"]["E1_mean"])
        e1_s = float(s2["DSA_Str_UQW1b_25"]["E1_mean"])
        ea_b = float(s2["DSA_Bent_UQW1b_25"]["E_amp_mean"])
        ea_s = float(s2["DSA_Str_UQW1b_25"]["E_amp_mean"])
        sm_b = float(s2["DSA_Bent_UQW1b_25"]["Smile_Frac_mean"])
        snr_s2 = float(s2["RAW_R0"]["SNR_img_dB_mean"])
        ea_red = 1.0 - ea_b / max(ea_s, 1e-12)
        x3_ok = (2.0 <= snr_s2 <= 6.0) and (((e1_b - e1_s) >= 0.005) or (ea_red >= 0.10)) and (sm_b <= 0.12)
        verdict["gates"]["X3_BentRay_Imaging_Benefit_S2"] = {
            "pass": bool(x3_ok),
            "S2_RAW_R0_SNR_dB": snr_s2,
            "E1_Bent_vs_Str": [e1_b, e1_s],
            "delta_E1": e1_b - e1_s,
            "Eamp_Bent_vs_Str": [ea_b, ea_s],
            "Eamp_red_pct": ea_red * 100.0,
            "Smile_Frac_Bent": sm_b,
        }
    else:
        verdict["gates"]["X3_BentRay_Imaging_Benefit_S2"] = {"pass": False, "status": "PENDING_WP2_S2"}

    # X4: Head-to-Head Uncertain-Q Showdown (S1-S4, overcompensation end m_lo)
    if len(s_data) == 4:
        snr_ok = all(2.0 <= float(s_data[sn]["W_deep"]["RAW_R0"]["SNR_img_dB_mean"]) <= 6.0 for sn in s_names)
        beat_std_cnt = 0
        beat_w1b_cnt = 0
        details = {}
        for sn in s_names:
            wd = s_data[sn]["W_deep"]
            p_lo = float(wd["DSA_Bent_UQW1b_25"]["E1"][0])
            std_max = max(float(wd["WE_QRTM_Unreg"]["E1"][0]), float(wd["WE_QRTM_ZhuLP"]["E1"][0]), float(wd["WE_QRTM_Reg"]["E1"][0]))
            we_w1b_lo = float(wd["DSA_WE_QRTM_W1b"]["E1"][0])
            ok_std = (p_lo > std_max)
            ok_w1b = (p_lo >= we_w1b_lo - 0.005)
            beat_std_cnt += int(ok_std)
            beat_w1b_cnt += int(ok_w1b)
            details[sn] = {
                "RAW_R0_SNR_dB": float(wd["RAW_R0"]["SNR_img_dB_mean"]),
                "m_lo_UQW1b": p_lo,
                "m_lo_max_std_WE": std_max,
                "m_lo_DSA_WE_W1b": we_w1b_lo,
                "delta_vs_max_std_WE": p_lo - std_max,
                "delta_vs_DSA_WE_W1b": p_lo - we_w1b_lo,
            }
        x4_ok = snr_ok and (beat_std_cnt == 4) and (beat_w1b_cnt >= 3)
        verdict["gates"]["X4_HardHeadToHead_UncertainQ"] = {
            "pass": bool(x4_ok),
            "all_snr_in_2_6_dB": snr_ok,
            "beat_all_std_WE_count": f"{beat_std_cnt}/4",
            "beat_or_tie_DSA_WE_W1b_count": f"{beat_w1b_cnt}/4",
            "details": details,
        }
    else:
        verdict["gates"]["X4_HardHeadToHead_UncertainQ"] = {"pass": False, "status": "PENDING_WP2"}

    # X5: Multipathing Crossover Audit on S3_SharpQ_Multipath
    if "S3_SharpQ_Multipath" in s_data:
        s3 = s_data["S3_SharpQ_Multipath"]
        nf3 = s3["NoiseFree_m1"]
        wd3 = s3["W_deep"]
        nf_we = float(nf3["WE_QRTM_Unreg"]["E1"])
        nf_uq = float(nf3["DSA_Bent_UQW1b_25"]["E1"])
        nf_hy = float(nf3["WE_Disp_UQW1b_25"]["E1"])
        min_uq = min(float(x) for x in wd3["DSA_Bent_UQW1b_25"]["E1"])
        min_hy = min(float(x) for x in wd3["WE_Disp_UQW1b_25"]["E1"])
        r0_s3  = float(wd3["RAW_R0"]["E1_mean"])
        x5_ok = (min_uq >= r0_s3 - 0.010) and (min_hy >= r0_s3 - 0.010)
        verdict["gates"]["X5_Multipath_Crossover_Audit_S3"] = {
            "pass": bool(x5_ok),
            "NF_m1_WE_QRTM_Unreg_vs_UQW1b_vs_Hybrid": [nf_we, nf_uq, nf_hy],
            "NF_m1_WE_advantage_over_ray_UQW1b": nf_we - nf_uq,
            "Noisy_min_E1_UQW1b_vs_Hybrid_vs_R0": [min_uq, min_hy, r0_s3],
            "Noisy_E1_mean_UQW1b_vs_Hybrid_vs_WE_Unreg_vs_WE_Reg_vs_DSA_WE_W1b": [
                float(wd3["DSA_Bent_UQW1b_25"]["E1_mean"]),
                float(wd3["WE_Disp_UQW1b_25"]["E1_mean"]),
                float(wd3["WE_QRTM_Unreg"]["E1_mean"]),
                float(wd3["WE_QRTM_Reg"]["E1_mean"]),
                float(wd3["DSA_WE_QRTM_W1b"]["E1_mean"]),
            ],
        }
    else:
        verdict["gates"]["X5_Multipath_Crossover_Audit_S3"] = {"pass": False, "status": "PENDING_WP2_S3"}

    # X6: Independent Time-Domain 8th-Order GSLS FD Verification (S5)
    s5_p = os.path.join(res_dir, "synth_S5_FD_GSLS_NoInvCrime.json")
    if os.path.isfile(s5_p):
        s5 = load_j(s5_p)["W_deep"]
        snr5 = float(s5["RAW_R0"]["SNR_img_dB_mean"])
        p_m = float(s5["DSA_Bent_UQW1b_25"]["E1_mean"])
        r0_m = float(s5["RAW_R0"]["E1_mean"])
        wu_m = float(s5["WE_QRTM_Unreg"]["E1_mean"])
        wr_m = float(s5["WE_QRTM_Reg"]["E1_mean"])
        p_lo = float(s5["DSA_Bent_UQW1b_25"]["E1"][0])
        ww_lo = float(s5["DSA_WE_QRTM_W1b"]["E1"][0])
        ea_red = 1.0 - float(s5["DSA_Bent_UQW1b_25"]["E_amp_mean"]) / max(float(s5["WE_Q_DispOnly"]["E_amp_mean"]), 1e-12)
        x6_ok = (2.0 <= snr5 <= 6.0) and (p_m > r0_m) and (p_m > wu_m) and (p_m > wr_m) and (p_lo > ww_lo) and (ea_red >= 0.25)
        verdict["gates"]["X6_Independent_GSLS_FD_S5"] = {
            "pass": bool(x6_ok),
            "S5_RAW_R0_SNR_dB": snr5,
            "E1_mean_UQW1b_vs_R0_WEUnreg_WEReg": [p_m, r0_m, wu_m, wr_m],
            "m05_UQW1b_vs_DSA_WE_W1b": [p_lo, ww_lo],
            "Eamp_red_vs_WE_DispOnly_pct": ea_red * 100.0,
        }
    else:
        verdict["gates"]["X6_Independent_GSLS_FD_S5"] = {"pass": False, "status": "PENDING_WP3"}

    # X7: Volve 4-25 Hz Real Field Showdown (WP4)
    v25_p = os.path.join(res_dir, "volve25_xmir_qrtm.json")
    if os.path.isfile(v25_p):
        v25 = load_j(v25_p)
        md = v25["methods"]
        prop = "DSA_Bent_UQW1b_PAnc_25"
        ghf = float(v25["G_HF_18_24"][prop])
        sm_m1 = float(md[prop]["Smile_Frac"][1])
        p_m05 = float(md[prop]["F5_R0band"][0])
        r0_f5 = float(md["RAW_R0"]["F5_R0band_mean"])
        we_u_m05 = float(md["WE_QRTM_Unreg_PAnc"]["F5_R0band"][0])
        x7_ok = (ghf >= 1.20) and (sm_m1 <= 0.08) and (p_m05 >= r0_f5 - 0.010)
        verdict["gates"]["X7_Volve25_Field_Showdown"] = {
            "pass": bool(x7_ok),
            "G_HF_18_24": ghf,
            "Smile_Frac_m1": sm_m1,
            "F5_Bent_UQW1b_mean": float(md[prop]["F5_R0band_mean"]),
            "F5_Str_UQW1b_mean": float(md["DSA_Str_UQW1b_PAnc_25"]["F5_R0band_mean"]),
            "F5_WE_Unreg_PAnc_mean": float(md["WE_QRTM_Unreg_PAnc"]["F5_R0band_mean"]),
            "F5_DSA_WE_Reg_PAnc_mean": float(md["DSA_WE_QRTM_Reg_PAnc"]["F5_R0band_mean"]),
            "F5_DSA_WE_W1b_PAnc_mean": float(md["DSA_WE_QRTM_W1b_PAnc"]["F5_R0band_mean"]),
            "m05_Bent_UQW1b_vs_WE_Unreg_vs_R0": [p_m05, we_u_m05, r0_f5],
        }
    else:
        verdict["gates"]["X7_Volve25_Field_Showdown"] = {"pass": False, "status": "PENDING_WP4"}

    n_pass = sum(1 for v in verdict["gates"].values() if v.get("pass"))
    n_tot = len(verdict["gates"])
    verdict["summary"] = {"passed": n_pass, "total": n_tot}
    out_p = os.path.join(res_dir, "xmir_qrtm_verdict.json")
    with open(out_p, "w", encoding="utf-8") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print(f"=== Pre-Registered Verdict (uq_rtm_xmir_qrtm_20261126): {n_pass}/{n_tot} PASS ===")
    for k, v in verdict["gates"].items():
        print(f"  {k:<40s}: {'PASS' if v.get('pass') else 'FAIL'}")
    print(f"Saved {out_p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
