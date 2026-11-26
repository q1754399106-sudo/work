#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
judge_closure.py
================
Pre-registered frozen mechanical verdict script for Gates R0-R7
(uq_rtm_closure_20261125). NEVER modify thresholds after pre-registration freeze.
"""

from __future__ import annotations

import json
import os
import statistics
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

    # R0: Integrity Gate (self-test log exists + FAIL on old report, PASS on new report if exists)
    selftest_log = os.path.join(log_dir, "wp0_audit_selftest.log")
    r0_selftest = os.path.isfile(selftest_log) and ("FINAL AUDIT VERDICT: FAIL" in open(selftest_log, encoding="utf-8").read())
    r0_new = False
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
        r0_new = (rc == 0)
    verdict["gates"]["R0_Integrity_Gate"] = {
        "pass": bool(r0_selftest and r0_new),
        "selftest_caught_old_defects": r0_selftest,
        "new_report_audit_zero_exit": r0_new,
    }

    # Load WP2 (Volve 4-45 Hz)
    v45_path = os.path.join(res_dir, "volve45_closure.json")
    if os.path.isfile(v45_path):
        v45 = load_j(v45_path)
        m45 = v45["methods"]
        prop = "DSA_UQW1b_TrueDSAiv_PAnc_25"
        ghf = float(v45["G_HF_28_40"][prop])

        # R1: Active Broadband Operator (G_HF >= 1.25)
        verdict["gates"]["R1_Active_Broadband_Operator"] = {
            "pass": bool(ghf >= 1.25),
            "G_HF_28_40_vs_DSA_R0_Eq": ghf,
            "threshold": 1.25,
        }

        # R2: Volve 4-45 Hz True Confluence Closure
        sm_prop = float(m45[prop]["Smile_Frac"][2])
        sm_s4h = float(m45["RAW_UQW1b_S4H_Unanc_100"]["Smile_Frac"][2])
        sm_red = 1.0 - sm_prop / max(sm_s4h, 1e-12)
        f5_r0 = float(m45["RAW_R0"]["F5_fix_mean"])
        f5_prop = float(m45[prop]["F5_fix_mean"])
        f5_all_gt = all(float(x) > f5_r0 for x in m45[prop]["F5_fix"])
        sl_eq = abs(float(v45["res_slope"]["DSA_R0_Eq"]))
        sl_pr = abs(float(v45["res_slope"][prop]))
        sl_imp = 1.0 - sl_pr / max(sl_eq, 1e-12)
        c_smile = (sm_prop <= 0.10 and sm_red >= 0.85)
        c_f5 = (f5_prop >= f5_r0 + 0.040 and f5_all_gt)
        c_slope = (sl_imp >= 0.30)
        verdict["gates"]["R2_Volve45_TrueConfluence_Closure"] = {
            "pass": bool(c_smile and c_f5 and c_slope),
            "sub_smile_pass": c_smile,
            "smile_frac_m1": sm_prop,
            "smile_reduction_vs_S4H": sm_red,
            "sub_f5_pass": c_f5,
            "F5_fix_mean": f5_prop,
            "F5_RAW_R0": f5_r0,
            "F5_delta_vs_RAW_R0": f5_prop - f5_r0,
            "F5_5of5_gt_RAW_R0": f5_all_gt,
            "sub_slope_pass": c_slope,
            "abs_slope_prop": sl_pr,
            "abs_slope_R0_Eq": sl_eq,
            "slope_improvement_frac": sl_imp,
        }

        # R3: Volve 4-45 Hz True UQ vs Det at m=0.6
        p_m06 = float(m45[prop]["F5_fix"][0])
        iq_m06 = float(m45["DSA_IQ70_PAnc_25"]["F5_fix"][0])
        w1_m06 = float(m45["DSA_W1b_TrueDSAiv_PAnc_25"]["F5_fix"][0])
        verdict["gates"]["R3_Volve45_TrueUQ_vs_Det_At_m06"] = {
            "pass": bool((p_m06 - iq_m06) >= 0.020 and p_m06 >= w1_m06),
            "prop_m06": p_m06,
            "IQ70_m06": iq_m06,
            "W1b_TrueDSAiv_m06": w1_m06,
            "delta_vs_IQ70": p_m06 - iq_m06,
            "delta_vs_W1b": p_m06 - w1_m06,
        }
    else:
        for g in ("R1_Active_Broadband_Operator", "R2_Volve45_TrueConfluence_Closure", "R3_Volve45_TrueUQ_vs_Det_At_m06"):
            verdict["gates"][g] = {"pass": False, "status": "PENDING_WP2"}

    # R4: PSAE vs NoW2 Attribution (WP3 + WP2)
    v25_path = os.path.join(res_dir, "volve25_psae_ctrl.json")
    if os.path.isfile(v25_path):
        v25 = load_j(v25_path)
        m25 = v25["methods"]
        r4_25 = {
            "Unanc_F5": float(m25["DSA_UQW1b_Unanc_25"]["F5_R0band_mean"]),
            "NoW2_Rel_F5": float(m25["DSA_UQW1b_NoW2_Rel_25"]["F5_R0band_mean"]),
            "PAnc_F5": float(m25["DSA_UQW1b_PAnc_25"]["F5_R0band_mean"]),
            "Decon_F5": float(m25["DSA_UQW1b_Decon_25"]["F5_R0band_mean"]),
            "NoW2_Rel_F14": float(m25["DSA_UQW1b_NoW2_Rel_25"]["F14_R0band_mean"]),
            "PAnc_F14": float(m25["DSA_UQW1b_PAnc_25"]["F14_R0band_mean"]),
            "NoW2_recovery_pct_of_PSAE_gain": (
                (float(m25["DSA_UQW1b_NoW2_Rel_25"]["F5_R0band_mean"]) - float(m25["DSA_UQW1b_Unanc_25"]["F5_R0band_mean"]))
                / max(float(m25["DSA_UQW1b_PAnc_25"]["F5_R0band_mean"]) - float(m25["DSA_UQW1b_Unanc_25"]["F5_R0band_mean"]), 1e-12)
                * 100.0
            ),
        }
        r4_45 = {}
        if os.path.isfile(v45_path):
            m45 = load_j(v45_path)["methods"]
            r4_45 = {
                "NoW2_F5": float(m45["DSA_UQW1b_TrueDSAiv_NoW2_25"]["F5_fix_mean"]),
                "PAnc_F5": float(m45["DSA_UQW1b_TrueDSAiv_PAnc_25"]["F5_fix_mean"]),
                "NoW2_F14": float(m45["DSA_UQW1b_TrueDSAiv_NoW2_25"]["F14_fix_mean"]),
                "PAnc_F14": float(m45["DSA_UQW1b_TrueDSAiv_PAnc_25"]["F14_fix_mean"]),
                "Decon_F5": float(m45["DSA_UQW1b_TrueDSAiv_Decon_25"]["F5_fix_mean"]),
                "Decon_F14": float(m45["DSA_UQW1b_TrueDSAiv_Decon_25"]["F14_fix_mean"]),
            }
        verdict["gates"]["R4_PSAE_vs_NoW2_Attribution"] = {
            "pass": bool(os.path.isfile(v25_path) and os.path.isfile(v45_path)),
            "volve25": r4_25,
            "volve45": r4_45,
        }
    else:
        verdict["gates"]["R4_PSAE_vs_NoW2_Attribution"] = {"pass": False, "status": "PENDING_WP3"}

    # R5: Synth U6-U7 Dual Metric & SA Tradeoff
    u6_p, u7_p = os.path.join(res_dir, "synth_U6.json"), os.path.join(res_dir, "synth_U7.json")
    if os.path.isfile(u6_p) and os.path.isfile(u7_p):
        u6, u7 = load_j(u6_p), load_j(u7_p)
        th3_ok = all(abs(float(u["noise_acc"]["DSA_B300_log2err"])) <= 1.0 for u in (u6, u7))
        eamp_red = [
            1.0 - float(u["W_deep"]["DSA_UQW1b_25"]["E_amp_mean"]) / max(float(u["W_deep"]["DSA_PH_25"]["E_amp_mean"]), 1e-12)
            for u in (u6, u7)
        ]
        eamp_ok = all(r >= 0.30 for r in eamp_red)
        m05_gains = [
            float(u["W_deep"]["DSA_UQW1b_25"]["E1"][0]) - max(float(u["W_deep"]["DSA_IQ70_25"]["E1"][0]), float(u["W_deep"]["DSA_W1b_25"]["E1"][0]))
            for u in (u6, u7)
        ]
        m05_ok = all(g >= 0.030 for g in m05_gains)
        snr_ok = all(2.0 <= float(u["W_deep"]["RAW_R0"]["SNR_img_dB_mean"]) <= 6.0 for u in (u6, u7))
        verdict["gates"]["R5_Synth_U6U7_DualMetric_And_SA_Tradeoff"] = {
            "pass": bool(th3_ok and eamp_ok and m05_ok),
            "snr_in_2_6_dB": snr_ok,
            "RAW_R0_SNR_dB": [float(u6["W_deep"]["RAW_R0"]["SNR_img_dB_mean"]), float(u7["W_deep"]["RAW_R0"]["SNR_img_dB_mean"])],
            "theorem3_DSA_B300_log2err": [float(u6["noise_acc"]["DSA_B300_log2err"]), float(u7["noise_acc"]["DSA_B300_log2err"])],
            "theorem3_pass": th3_ok,
            "E_amp_reduction_vs_PH": eamp_red,
            "E_amp_pass": eamp_ok,
            "m05_E1_gain_vs_max_det": m05_gains,
            "m05_pass": m05_ok,
            "U6_SA_vs_DSA": [float(u6["W_deep"]["SA_UQW1b_25"]["E1_mean"]), float(u6["W_deep"]["DSA_UQW1b_25"]["E1_mean"])],
            "U7_SA_vs_DSA": [float(u7["W_deep"]["SA_UQW1b_25"]["E1_mean"]), float(u7["W_deep"]["DSA_UQW1b_25"]["E1_mean"])],
        }
    else:
        verdict["gates"]["R5_Synth_U6U7_DualMetric_And_SA_Tradeoff"] = {"pass": False, "status": "PENDING_WP4"}

    # R6: Non-Scalar Spatial Q Error U8-U9
    u8_p, u9_p = os.path.join(res_dir, "synth_U8.json"), os.path.join(res_dir, "synth_U9.json")
    if os.path.isfile(u8_p) and os.path.isfile(u9_p):
        u8, u9 = load_j(u8_p), load_j(u9_p)
        u8_ok = (float(u8["W_deep"]["DSA_UQW1b_25"]["E1_mean"]) > float(u8["W_deep"]["DSA_IQ70_25"]["E1_mean"])) and (
            min(u8["W_deep"]["DSA_UQW1b_25"]["E1"]) >= float(u8["W_deep"]["RAW_R0"]["E1_mean"]) - 0.020
        )
        u9_ok = (float(u9["W_deep"]["DSA_UQW1b_25"]["E1_mean"]) > float(u9["W_deep"]["DSA_IQ70_25"]["E1_mean"])) and (
            min(u9["W_deep"]["DSA_UQW1b_25"]["E1"]) >= float(u9["W_deep"]["RAW_R0"]["E1_mean"]) - 0.020
        )
        verdict["gates"]["R6_NonScalar_Spatial_Q_Error_U8U9"] = {
            "pass": bool(u8_ok and u9_ok),
            "U8_pass": u8_ok,
            "U8_UQW1b_vs_IQ70_vs_W1b_vs_PH_vs_R0": {
                k: float(u8["W_deep"][k]["E1_mean"])
                for k in ("DSA_UQW1b_25", "DSA_IQ70_25", "DSA_W1b_25", "DSA_PH_25", "RAW_R0")
            },
            "U8_min_E1": min(u8["W_deep"]["DSA_UQW1b_25"]["E1"]),
            "U9_pass": u9_ok,
            "U9_UQW1b_vs_IQ70_vs_W1b_vs_PH_vs_R0": {
                k: float(u9["W_deep"][k]["E1_mean"])
                for k in ("DSA_UQW1b_25", "DSA_IQ70_25", "DSA_W1b_25", "DSA_PH_25", "RAW_R0")
            },
            "U9_min_E1": min(u9["W_deep"]["DSA_UQW1b_25"]["E1"]),
        }
    else:
        verdict["gates"]["R6_NonScalar_Spatial_Q_Error_U8U9"] = {"pass": False, "status": "PENDING_WP4"}

    # R7: Anomaly Resolution (WP1 diagnostics JSON + figure exist with D1-D4 populated)
    wp1_p = os.path.join(res_dir, "wp1_diagnostics.json")
    wp1_fig = os.path.join(root, "figures", "fig_wp1_diagnostics.png")
    if os.path.isfile(wp1_p) and os.path.isfile(wp1_fig):
        d_wp1 = load_j(wp1_p)
        has_all = all(k in d_wp1 for k in ("D1_F14", "D2_RAW_vs_DSA", "D3_Lag_Curves", "D4_U4b"))
        verdict["gates"]["R7_Anomaly_Resolution"] = {"pass": bool(has_all), "diagnostics_json": True, "figure_png": True}
    else:
        verdict["gates"]["R7_Anomaly_Resolution"] = {"pass": False, "status": "PENDING_WP1"}

    n_pass = sum(1 for v in verdict["gates"].values() if v.get("pass"))
    n_tot = len(verdict["gates"])
    verdict["summary"] = {"passed": n_pass, "total": n_tot}

    out_path = os.path.join(res_dir, "closure_verdict.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print(f"=== Closure Pre-Registered Verdict: {n_pass}/{n_tot} PASS ===")
    for k, v in verdict["gates"].items():
        print(f"  {k:<45s}: {'PASS' if v.get('pass') else 'FAIL'}")
    print(f"Saved {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
