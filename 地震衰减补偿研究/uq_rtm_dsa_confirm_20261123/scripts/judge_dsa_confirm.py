"""Pre-registered mechanical confirmation judge for DSA-UQW1 (uq_rtm_dsa_confirm_20261123/PLAN_ZH.md)."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
MODELS = ("H1", "H2", "H3", "H4")


def load_j(name):
    p = RES / f"{name}.json"
    if not p.exists():
        raise SystemExit(f"Missing {p}")
    return json.loads(p.read_text(encoding="utf-8"))


def vec(x):
    return np.atleast_1d(np.asarray(x, dtype=float))


def main():
    D = {m: load_j(m) for m in MODELS}
    H5 = load_j("H5_VolveAnc")

    # K1_DSA_vs_RAW on noisy W_deep
    dE1, dSNR, dHF, dAmp = {}, {}, {}, {}
    sa_dE1, sa_dSNR = {}, {}
    for m in MODELS:
        wd = D[m]["noisy"]["W_deep"]
        dE1[m]  = float(wd["DSA_UQW1"]["E1_mean"] - wd["RAW_UQW1"]["E1_mean"])
        dSNR[m] = float(wd["DSA_UQW1"]["SNR_img_dB_mean"] - wd["RAW_UQW1"]["SNR_img_dB_mean"])
        dHF[m]  = float(wd["DSA_UQW1"]["E1_HF_mean"] - wd["RAW_UQW1"]["E1_HF_mean"])
        dAmp[m] = float(wd["DSA_UQW1"]["E_amp_mean"] - wd["RAW_UQW1"]["E_amp_mean"])
        sa_dE1[m]  = float(wd["SA_flat_UQW1"]["E1_mean"] - wd["RAW_UQW1"]["E1_mean"])
        sa_dSNR[m] = float(wd["SA_flat_UQW1"]["SNR_img_dB_mean"] - wd["RAW_UQW1"]["SNR_img_dB_mean"])

    k1_ok = bool(
        np.mean(list(dE1.values())) >= 0.020
        and all(v >= 0.010 for v in dE1.values())
        and np.mean(list(dSNR.values())) >= 1.00
        and all(v >= 0.0 for v in dSNR.values())
    )

    # K2_DSA_Fixes_SteepDip on H2
    h2_n = D["H2"]["noisy"]["W_deep"]
    h2_c = D["H2"]["clean"]["W_deep"]
    k2_det = dict(
        noisy_DSA_minus_SAflat_E1=float(h2_n["DSA_UQW1"]["E1_mean"] - h2_n["SA_flat_UQW1"]["E1_mean"]),
        noisy_DSA_minus_SAflat_HF=float(h2_n["DSA_UQW1"]["E1_HF_mean"] - h2_n["SA_flat_UQW1"]["E1_HF_mean"]),
        noisy_DSA_minus_SAflat_SNR=float(h2_n["DSA_UQW1"]["SNR_img_dB_mean"] - h2_n["SA_flat_UQW1"]["SNR_img_dB_mean"]),
        clean_DSA_minus_RAW_E1=float(h2_c["DSA_UQW1"]["E1_mean"] - h2_c["RAW_UQW1"]["E1_mean"]),
        clean_DSA_minus_RAW_HF=float(h2_c["DSA_UQW1"]["E1_HF_mean"] - h2_c["RAW_UQW1"]["E1_HF_mean"]),
        noisy_DSA_minus_RAW_E1=float(h2_n["DSA_UQW1"]["E1_mean"] - h2_n["RAW_UQW1"]["E1_mean"]),
        noisy_DSA_minus_RAW_HF=float(h2_n["DSA_UQW1"]["E1_HF_mean"] - h2_n["RAW_UQW1"]["E1_HF_mean"]),
    )
    k2_ok = bool(
        k2_det["noisy_DSA_minus_SAflat_E1"] >= 0.008
        and k2_det["noisy_DSA_minus_SAflat_HF"] >= 0.010
        and k2_det["noisy_DSA_minus_SAflat_SNR"] >= 0.50
        and k2_det["clean_DSA_minus_RAW_E1"] >= 0.0
        and k2_det["noisy_DSA_minus_RAW_E1"] >= 0.0
        and k2_det["clean_DSA_minus_RAW_HF"] >= -0.008
        and k2_det["noisy_DSA_minus_RAW_HF"] >= -0.008
    )

    # K3_Fair_UQ_Robustness on H1-H4 noisy W_deep
    k3_tbl = {}
    for m in MODELS:
        wd = D[m]["noisy"]["W_deep"]
        ms = vec(D[m]["m"])
        r0 = np.full(len(ms), float(vec(wd["DSA_R0"]["E1"])[0]))
        iq = vec(wd["DSA_IQ70"]["E1"])
        w1 = vec(wd["DSA_W1"]["E1"])
        uq = vec(wd["DSA_UQW1"]["E1"])
        bplus = np.maximum.reduce([r0, iq, w1])
        d = uq - bplus
        m05_vs_comp = float(uq[0] - max(iq[0], w1[0]))
        M2 = float(d[ms == 1.0][0])
        M3 = float(np.min(d[ms > 1.0]))
        wc_beats_comp = bool(np.min(uq) > np.min(iq) and np.min(uq) > np.min(w1))
        k3_tbl[m] = dict(
            m05_vs_comp=m05_vs_comp, M2=M2, M3=M3,
            min_uq=float(np.min(uq)), min_w1=float(np.min(w1)), min_iq=float(np.min(iq)),
            worst_beats_comp=wc_beats_comp,
        )
    k3_ok = bool(
        np.mean([k3_tbl[m]["m05_vs_comp"] for m in MODELS]) >= 0.050
        and all(k3_tbl[m]["m05_vs_comp"] > 0.0 for m in MODELS)
        and all(k3_tbl[m]["M2"] >= -0.025 for m in MODELS)
        and all(k3_tbl[m]["M3"] >= -0.035 for m in MODELS)
        and all(k3_tbl[m]["worst_beats_comp"] for m in MODELS)
    )

    # K4_Amp_NoHarm on H1-H4 noisy W_deep
    k4_ok = bool(all(v <= 0.020 for v in dAmp.values()) and np.mean(list(dAmp.values())) <= 0.0)

    # K5_VolveAnc_WellTie on H5_VolveAnc
    anc_t = vec(H5["DSA_UQW1_Anc"]["tie"])
    raw_t = vec(H5["RAW_UQW1"]["tie"])
    r0_t = float(H5["RAW_R0"])
    k5_det = dict(
        RAW_R0=r0_t,
        DSA_R0=float(H5["DSA_R0"]),
        RAW_UQW1_mean=float(np.mean(raw_t)),
        RAW_UQW1_m1=float(raw_t[2]),
        DSA_UQW1_Anc_mean=float(np.mean(anc_t)),
        DSA_UQW1_Anc_m1=float(anc_t[2]),
        d_Anc_minus_RAW_UQ_mean=float(np.mean(anc_t) - np.mean(raw_t)),
        d_Anc_minus_RAW_UQ_m1=float(anc_t[2] - raw_t[2]),
        d_Anc_minus_RAW_R0_mean=float(np.mean(anc_t) - r0_t),
        d_Anc_minus_RAW_R0_m1=float(anc_t[2] - r0_t),
    )
    k5_ok = bool(
        k5_det["d_Anc_minus_RAW_UQ_mean"] >= 0.050
        and k5_det["d_Anc_minus_RAW_UQ_m1"] >= 0.050
        and k5_det["d_Anc_minus_RAW_R0_mean"] >= 0.0
        and k5_det["d_Anc_minus_RAW_R0_m1"] >= 0.0
    )

    verdict = dict(
        K1_DSA_vs_RAW=dict(
            pass_rule=k1_ok,
            mean_dE1=float(np.mean(list(dE1.values()))),
            mean_dSNR_dB=float(np.mean(list(dSNR.values()))),
            by_model_dE1=dE1, by_model_dSNR=dSNR, by_model_dHF=dHF,
            sa_flat_dE1=sa_dE1, sa_flat_dSNR=sa_dSNR,
        ),
        K2_DSA_Fixes_SteepDip=dict(pass_rule=k2_ok, H2=k2_det),
        K3_Fair_UQ_Robustness=dict(pass_rule=k3_ok, by_model=k3_tbl),
        K4_Amp_NoHarm=dict(pass_rule=k4_ok, mean_dAmp=float(np.mean(list(dAmp.values()))), by_model_dAmp=dAmp),
        CORE_CONFIRM_PASS=bool(k1_ok and k2_ok and k3_ok and k4_ok),
        K5_VolveAnc_WellTie=dict(pass_rule=k5_ok, details=k5_det),
    )
    out = RES / "dsa_confirm_verdict.json"
    out.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
