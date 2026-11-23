"""Pre-registered mechanical judge for SPG-UQ-RTM (uq_rtm_spg_20261123/PLAN_ZH.md)."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
MODELS = ("G1", "G2", "G3", "G4")


def load_model(m):
    p = RES / f"{m}.json"
    if not p.exists():
        raise SystemExit(f"Missing {p}")
    return json.loads(p.read_text(encoding="utf-8"))


def vec(x):
    return np.atleast_1d(np.asarray(x, dtype=float))


def main():
    D = {m: load_model(m) for m in MODELS}

    # Gate 1: G1_SPG_vs_RAW on noisy W_deep
    dE1_spg_raw = {}
    dSNR_spg_raw = {}
    dE1_sa_raw = {}
    dSNR_sa_raw = {}
    for m in MODELS:
        wd = D[m]["noisy"]["W_deep"]
        dE1_spg_raw[m] = float(wd["SPG_UQW1"]["E1_mean"] - wd["RAW_UQW1"]["E1_mean"])
        dSNR_spg_raw[m] = float(wd["SPG_UQW1"]["SNR_img_dB_mean"] - wd["RAW_UQW1"]["SNR_img_dB_mean"])
        dE1_sa_raw[m] = float(wd["SA_UQW1"]["E1_mean"] - wd["RAW_UQW1"]["E1_mean"])
        dSNR_sa_raw[m] = float(wd["SA_UQW1"]["SNR_img_dB_mean"] - wd["RAW_UQW1"]["SNR_img_dB_mean"])

    e1_vals = list(dE1_spg_raw.values())
    snr_vals = list(dSNR_spg_raw.values())
    g1_ok = bool(
        np.mean(e1_vals) >= 0.025
        and sum(v >= 0.015 for v in e1_vals) >= 3
        and all(v >= 0.0 for v in e1_vals)
        and np.mean(snr_vals) >= 0.80
        and all(v > 0.0 for v in snr_vals)
    )

    # Gate 2: G2_Fair_UQ_Margin (SPG_UQW1 vs SPG_B+ = max(SPG_R0, SPG_IQ70, SPG_W1)) on noisy W_deep
    fair_tbl = {}
    for m in MODELS:
        wd = D[m]["noisy"]["W_deep"]
        ms = vec(D[m]["m"])
        r0 = np.full(len(ms), float(vec(wd["SPG_R0"]["E1"])[0]))
        iq = vec(wd["SPG_IQ70"]["E1"])
        w1 = vec(wd["SPG_W1"]["E1"])
        uq = vec(wd["SPG_UQW1"]["E1"])
        bplus = np.maximum.reduce([r0, iq, w1])
        d = uq - bplus
        M1 = float(np.mean(d[ms < 1.0]))
        M2 = float(d[ms == 1.0][0])
        M3 = float(np.min(d[ms > 1.0]))
        # Also compare RAW_UQW1 - RAW_W1 vs SPG_UQW1 - SPG_W1 at m=0.5
        raw_m05_margin = float(vec(wd["RAW_UQW1"]["E1"])[0] - vec(wd["RAW_W1"]["E1"])[0])
        spg_m05_margin = float(uq[0] - w1[0])
        fair_tbl[m] = dict(M1=M1, M2=M2, M3=M3, raw_m05_uq_minus_w1=raw_m05_margin, spg_m05_uq_minus_w1=spg_m05_margin)

    g2_ok = bool(
        np.mean([fair_tbl[m]["M1"] for m in MODELS]) >= 0.050
        and all(fair_tbl[m]["M1"] >= 0.030 for m in MODELS)
        and all(fair_tbl[m]["M2"] >= -0.025 for m in MODELS)
        and all(fair_tbl[m]["M3"] >= -0.040 for m in MODELS)
    )

    # Gate 3: G3_SteepDip_Safety on G2 (clean & noisy W_deep)
    g3_det = {}
    for rg in ("clean", "noisy"):
        wd = D["G2"][rg]["W_deep"]
        g3_det[rg] = dict(
            dE1=float(wd["SPG_UQW1"]["E1_mean"] - wd["RAW_UQW1"]["E1_mean"]),
            dE1_HF=float(wd["SPG_UQW1"]["E1_HF_mean"] - wd["RAW_UQW1"]["E1_HF_mean"]),
            dSNR=float(wd["SPG_UQW1"]["SNR_img_dB_mean"] - wd["RAW_UQW1"]["SNR_img_dB_mean"]),
        )
    g3_ok = bool(all(g3_det[rg]["dE1"] >= -0.005 and g3_det[rg]["dE1_HF"] >= -0.005 for rg in ("clean", "noisy")))

    # Gate 4: G4_Chimney_Amp_Safety on G1 & G4 noisy W_chim
    g4_det = {}
    for m in ("G1", "G4"):
        wc = D[m]["noisy"]["W_chim"]
        g4_det[m] = dict(
            raw_E_amp=float(wc["RAW_UQW1"]["E_amp_mean"]),
            sa_E_amp=float(wc["SA_UQW1"]["E_amp_mean"]),
            spg_E_amp=float(wc["SPG_UQW1"]["E_amp_mean"]),
            d_spg_minus_raw=float(wc["SPG_UQW1"]["E_amp_mean"] - wc["RAW_UQW1"]["E_amp_mean"]),
        )
    g4_ok = bool(all(g4_det[m]["d_spg_minus_raw"] <= 0.030 for m in ("G1", "G4")))

    verdict = dict(
        G1_SPG_vs_RAW=dict(
            pass_rule=g1_ok,
            mean_dE1=float(np.mean(e1_vals)),
            mean_dSNR_dB=float(np.mean(snr_vals)),
            by_model_dE1=dE1_spg_raw,
            by_model_dSNR=dSNR_spg_raw,
            sa_only_dE1=dE1_sa_raw,
            sa_only_dSNR=dSNR_sa_raw,
        ),
        G2_Fair_UQ_Margin=dict(pass_rule=g2_ok, by_model=fair_tbl),
        G3_SteepDip_Safety=dict(pass_rule=g3_ok, G2=g3_det),
        G4_Chimney_Amp_Safety=dict(pass_rule=g4_ok, by_model=g4_det),
        ALL_PASS=bool(g1_ok and g2_ok and g3_ok and g4_ok),
    )
    out = RES / "spg_verdict.json"
    out.write_text(json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
