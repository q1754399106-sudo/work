"""Shared code of experiment 16 (PLAN_ZH.md): does dropping SC-a's lateral coherence term (mu_relative 10 -> 0)
help on unseen faulted, complex-Q data, for the SR route, the user's QScan route and the true effective Q?

Imported unchanged: lib15 (experiment 15; through it the experiment-7 to -13 code) and band_last of the SC-a
ceiling check (descriptive only). New code: data loading of the second batch, the method set.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "fault_q_confirm_20261030" / "scripts"))
sys.path.insert(0, str(BASE / "sca_ceiling_20261030" / "scripts"))
import lib15  # noqa: E402
from lib15 import np, lib, semi, sq, sc11, lib12, TW, SEL, A_TUNED, G_IQ, score, masks, local  # noqa: E402,F401
from libc import band_last  # noqa: E402

DATA = BASE / "semisynth_fault_q2_20261031"
META = json.loads((DATA / "analysis" / "datasets.json").read_text(encoding="utf-8"))
TEST = [f"Z{k}_c{c}_r{r}" for k in "DEF" for c in (0.75, 1.05) for r in (0, 1)]
SMOKE = ["ZA_c0.75_r0"]                                   # development structure A of the first batch
A_ORIG, MU = lib.SEL_SC                                    # 0.125, 10 (the user's SC-a)


def load(did):
    if did not in META:
        return lib15.load(did)
    return {k: v.astype(float) for k, v in np.load(DATA / "cache" / "data" / f"{did}.npz").items()}, META[did]


def methods(D, cm, Q, c_sr, c_qs, c_or, p):
    """PLAN_ZH.md section 3. c_sr: SR factor; c_qs: QScan factor; c_or: c_eff_cum of the truth."""
    sc = lambda c, a, mu: lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], a, mu)[0]
    out = dict(Input=D, IQ70=lib12.iq(D, cm, Q, c_sr, c_sr, G_IQ))
    out["VDB"] = lib12.tvw(lib12.iq(D, cm, Q, None, c_sr, 0), p["wv"]["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"])
    out["SR_mu10"] = sc(c_sr, A_TUNED, MU); out["SR_mu1"] = sc(c_sr, A_TUNED, 1.0); out["SR_mu0"] = sc(c_sr, A_TUNED, 0.0)
    out["QSCa"] = sc(c_qs, A_ORIG, MU); out["QSCa_mu0"] = sc(c_qs, A_ORIG, 0.0)
    out["QT_mu10"] = sc(c_qs, A_TUNED, MU); out["QT_mu0"] = sc(c_qs, A_TUNED, 0.0)
    out["OR_mu10"] = sc(c_or, A_TUNED, MU); out["OR_mu0"] = sc(c_or, A_TUNED, 0.0)
    return out
