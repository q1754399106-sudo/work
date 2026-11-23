"""Pre-registered mechanical judge for T2.1 (wsuq_rtm_20261123/PLAN_T21_ZH.md)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EIK = ROOT / "results" / "wsuq_eikonal_bench.json"
RTM = ROOT / "results" / "wsuq_rtm2d_bench.json"


def main():
    eik = json.loads(EIK.read_text(encoding="utf-8"))
    rtm = json.loads(RTM.read_text(encoding="utf-8"))
    verdict = {}

    # 1. T21_Eik_Volve (Volve F-5 124-shot exact TS_abs, 4-45 Hz, 3 dB)
    v3 = eik["Volve_F5_SNR3dB"]
    d_snr_mean = v3["WS_UQ1"]["SNR_dB_mean"] - v3["UQW1_old"]["SNR_dB_mean"]
    d_snr_m05 = v3["WS_UQ1"]["SNR_dB"][0] - v3["UQW1_old"]["SNR_dB"][0]
    d_e1hf = v3["WS_UQ1"]["E1_HF_mean"] - v3["UQW1_old"]["E1_HF_mean"]
    t21_eik_ok = (d_snr_mean >= 1.0) and (d_snr_m05 >= 2.5) and (d_e1hf > 0.0)
    verdict["T21_Eik_Volve"] = dict(
        pass_rule=bool(t21_eik_ok),
        d_SNR_mean_dB=d_snr_mean, d_SNR_m05_dB=d_snr_m05,
        d_E1_HF_mean=d_e1hf,
        WS_UQ1_SNR_mean=v3["WS_UQ1"]["SNR_dB_mean"], UQW1_old_SNR_mean=v3["UQW1_old"]["SNR_dB_mean"],
        WS_VA_UQ1_SNR_mean=v3["WS_VA_UQ1"]["SNR_dB_mean"], VA_old_SNR_mean=v3["VA_old"]["SNR_dB_mean"],
    )

    # 2. T21_RTM2D_Gain (true 2D full-wave gas-chimney RTM, noisy 3 dB, W_chim)
    nc = rtm["noisy"]["W_chim"]
    d_hf_rtm = nc["WS_UQ1"]["E1_HF_mean"] - nc["UQW1_old"]["E1_HF_mean"]
    d_snr_rtm = nc["WS_UQ1"]["SNR_img_dB_mean"] - nc["UQW1_old"]["SNR_img_dB_mean"]
    d_m05_rtm = nc["WS_UQ1"]["E1"][0] - nc["UQW1_old"]["E1"][0]
    t21_gain_ok = (d_hf_rtm > 0.0) and (d_snr_rtm >= 0.5) and (d_m05_rtm > 0.0)
    verdict["T21_RTM2D_Gain"] = dict(
        pass_rule=bool(t21_gain_ok),
        d_E1_HF_mean=d_hf_rtm, d_SNR_img_dB_mean=d_snr_rtm, d_E1_m05=d_m05_rtm,
        WS_UQ1_E1_mean=nc["WS_UQ1"]["E1_mean"], UQW1_old_E1_mean=nc["UQW1_old"]["E1_mean"],
        WS_UQ1_E1_HF_mean=nc["WS_UQ1"]["E1_HF_mean"], UQW1_old_E1_HF_mean=nc["UQW1_old"]["E1_HF_mean"],
        WS_UQ1_SNR_dB_mean=nc["WS_UQ1"]["SNR_img_dB_mean"], UQW1_old_SNR_dB_mean=nc["UQW1_old"]["SNR_img_dB_mean"],
        WS_VA_UQ1_SNR_dB_mean=nc["WS_VA_UQ1"]["SNR_img_dB_mean"], VA_old_SNR_dB_mean=nc["VA_old"]["SNR_img_dB_mean"],
    )

    # 3. T21_RTM2D_NoHarm (clean 15 dB, m=1.0 index 2)
    cc = rtm["clean"]["W_chim"]
    d_clean_m1 = cc["WS_UQ1"]["E1"][2] - cc["UQW1_old"]["E1"][2]
    verdict["T21_RTM2D_NoHarm"] = dict(
        pass_rule=bool(d_clean_m1 >= -0.01),
        d_E1_clean_m1=d_clean_m1,
        WS_UQ1_E1_m1=cc["WS_UQ1"]["E1"][2], UQW1_old_E1_m1=cc["UQW1_old"]["E1"][2],
    )

    out_path = ROOT / "results" / "t21_verdict.json"
    out_path.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
