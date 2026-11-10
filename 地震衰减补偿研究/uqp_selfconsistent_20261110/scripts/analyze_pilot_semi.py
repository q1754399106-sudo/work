"""Tables of pilot_semi.py (dataset means; 90% bootstrap intervals over datasets, 2000 draws, seed 20261125).
Descriptive. usage: analyze_pilot_semi.py"""
import json, glob
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; rng = np.random.default_rng(20261125)
KEYS = ("IQ70", "W0", "UQW", "MW0_K1", "UQMW_K1", "MW0_K5", "UQMW_K5", "MW0_K9", "UQMW_K9")


def ci(d):
    b = rng.integers(0, len(d), (2000, len(d))); m = d[b].mean(axis=1); return d.mean(), np.quantile(m, 0.05), np.quantile(m, 0.95), (d > 0).sum(), len(d)


for v in ("cal", "nd"):
    R = [json.load(open(p, encoding="utf-8")) for p in sorted(glob.glob(str(ROOT / "results" / "pilot_semi" / f"*_{v}.json")))]
    if not R:
        continue
    ms = [r["m"] for r in R[0]["records"]]; i1 = ms.index(1.0); T = {k: np.array([[r[k] for r in d["records"]] for d in R]) for k in KEYS}
    print(f"\n===== {v}: {len(R)} datasets\n method     " + "  ".join("%5.2f" % m for m in ms) + "    mean")
    for k in KEYS:
        x = T[k].mean(axis=0); print(" %-9s  " % k + "  ".join("%.3f" % a for a in x) + "   %.3f" % x.mean())
    sr = {k: np.array([d["sr"][k] for d in R]) for k in ("IQ70", "W0", "MW0_K1", "MW0_K5", "MW0_K9", "ratio")}
    print(" SR route: c_sr / c_true median %.3f, |log| median %.3f, max %.3f" % (np.median(sr["ratio"]), np.median(np.abs(np.log(sr["ratio"]))), np.abs(np.log(sr["ratio"])).max()))
    print(" SR route scores: " + "  ".join("%s %.3f" % (k, sr[k].mean()) for k in ("IQ70", "W0", "MW0_K1", "MW0_K5", "MW0_K9")))
    for name, d in (("scenario mean UQMW_K9 - UQW", T["UQMW_K9"].mean(axis=1) - T["UQW"].mean(axis=1)),
                    ("scenario mean UQMW_K5 - UQW", T["UQMW_K5"].mean(axis=1) - T["UQW"].mean(axis=1)),
                    ("scenario mean MW0_K9 - W0", T["MW0_K9"].mean(axis=1) - T["W0"].mean(axis=1)),
                    ("scenario mean MW0_K9 - IQ70", T["MW0_K9"].mean(axis=1) - T["IQ70"].mean(axis=1)),
                    ("scenario mean UQMW_K9 - MW0_K9", T["UQMW_K9"].mean(axis=1) - T["MW0_K9"].mean(axis=1)),
                    ("m = 1: MW0_K9 - W0", T["MW0_K9"][:, i1] - T["W0"][:, i1]),
                    ("m = 1: MW0_K9 - IQ70", T["MW0_K9"][:, i1] - T["IQ70"][:, i1]),
                    ("m = 1: MW0_K5 - MW0_K1", T["MW0_K5"][:, i1] - T["MW0_K1"][:, i1]),
                    ("SR route: MW0_K9 - IQ70", sr["MW0_K9"] - sr["IQ70"]), ("SR route: MW0_K9 - W0", sr["MW0_K9"] - sr["W0"]),
                    ("SR route: MW0_K5 - IQ70", sr["MW0_K5"] - sr["IQ70"]),
                    ("SR-route MW0_K9 - scenario-mean UQMW_K9", sr["MW0_K9"] - T["UQMW_K9"].mean(axis=1)),
                    ("steering: MW0_K9 - no steering (m = 1)", T["MW0_K9"][:, i1] - np.array([d["diag"]["MW0_K9_nosteer"] for d in R]))):
        print(" %-42s %+.3f [%+.3f, %+.3f]  positive %d/%d" % ((name,) + ci(d)))
    if v == "cal":                                    # archived SC-a (a 2^-13, mu 0.25) on the same data: stage 0 and experiment 17
        Q0 = ROOT.parent / "qunknown_comp_20261109" / "results" / "stage0"; E17 = ROOT.parent / "mu_tune_20261031" / "results" / "dev"
        msc = (0.5, 0.7, 1.0, 1.4, 2.0); js = [ms.index(m) for m in msc]
        sca = np.array([[json.load(open(Q0 / f"{d['id']}__SCa_m{m:g}.json", encoding="utf-8"))["score"]["S"] for m in msc] for d in R])
        scsr = np.array([json.load(open(E17 / f"{d['id']}.json", encoding="utf-8"))["scores"]["SR_mu0.25"]["S"] for d in R])
        print(" archived SC-a by m " + "  ".join("%g: %.3f" % (m, x) for m, x in zip(msc, sca.mean(axis=0))) + "   mean %.3f;  SR route %.3f" % (sca.mean(), scsr.mean()))
        for name, d in (("5 m mean: UQMW_K9 - SC-a", T["UQMW_K9"][:, js].mean(axis=1) - sca.mean(axis=1)),
                        ("5 m mean: MW0_K9 - SC-a", T["MW0_K9"][:, js].mean(axis=1) - sca.mean(axis=1)),
                        ("m = 1: MW0_K9 - SC-a", T["MW0_K9"][:, i1] - sca[:, 2]), ("SR route: MW0_K9 - SC-a", sr["MW0_K9"] - scsr)):
            print(" %-42s %+.3f [%+.3f, %+.3f]  positive %d/%d" % ((name,) + ci(d)))
    print(" coherence q (10-60 Hz median): steered %.3f, not steered %.3f; |dip| median %.2f ms/trace, p95 %.2f" % (
        np.mean([d["diag"]["q_steer"] for d in R]), np.mean([d["diag"]["q_nosteer"] for d in R]),
        np.mean([d["diag"]["dip_abs_med_ms"] for d in R]), np.mean([d["diag"]["dip_abs_p95_ms"] for d in R])))
