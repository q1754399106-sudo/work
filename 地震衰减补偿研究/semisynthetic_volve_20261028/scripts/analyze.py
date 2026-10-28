"""Experiment 8 summary and hypotheses H1-H5 (PLAN_ZH.md) from results/scores/*.json."""
import sys, json, csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import semi  # noqa: E402
from semi import np, bm, C_TRUE  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = semi.ROOT; SC = ROOT / "results" / "scores"; AN = ROOT / "analysis"
METHODS = ["Input", "IQ10", "IQ20", "IQ30", "IQ20a", "QIQa", "SC", "SCa", "QSCa", "SCa_oracle"]
KEYS = ["M1", "M1w"] + [f"M2_{a}_{b}" for a, b in semi.BANDS] + ["T2like", "band95_Hz"]


def load():
    return [json.loads(p.read_text(encoding="utf-8")) for p in sorted(SC.glob("*.json"))]


def ci(d):
    d = np.asarray(d, float); m = d.mean()
    if len(d) < 2:
        return [float(m), float("nan"), float("nan")]
    h = stats.t.ppf(0.95, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))
    return [float(m), float(m - h), float(m + h)]


def group(rows, g, c):
    return [r for r in rows if r["id"].startswith(g + "_") and abs(r["c_true"] - c) < 1e-9]


def diff(rs, a, b, key="M1"):
    return [r["scores"][a][key] - r["scores"][b][key] for r in rs]


def main():
    rows = load(); out = {"n_datasets": len(rows)}; table = []
    for g in ("P", "G", "R1", "PH"):
        for c in C_TRUE:
            rs = group(rows, g, c)
            if not rs:
                continue
            for m in METHODS:
                table.append(dict(group=g, c_true=c, method=m, n=len(rs),
                                  **{k: float(np.mean([r["scores"][m][k] for r in rs])) for k in KEYS}))
    with open(AN / "table.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0])); w.writeheader(); w.writerows(table)
    H = {}
    for c in C_TRUE:
        rs = group(rows, "P", c); h = {"n": len(rs)}
        h["H1_QSCa_minus_Input"] = ci(diff(rs, "QSCa", "Input")); h["H1_pass"] = h["H1_QSCa_minus_Input"][1] > 0
        h["SCa_minus_Input"] = ci(diff(rs, "SCa", "Input"))
        h["H2_QSCa_minus_IQ20a"] = ci(diff(rs, "QSCa", "IQ20a")); h["H2_QSCa_minus_QIQa"] = ci(diff(rs, "QSCa", "QIQa"))
        h["QSCa_minus_SC"] = ci(diff(rs, "QSCa", "SC")); h["QSCa_minus_oracle"] = ci(diff(rs, "QSCa", "SCa_oracle"))
        bands = {}
        for a, b in semi.BANDS:
            k = f"M2_{a}_{b}"; mi = float(np.mean([r["scores"]["Input"][k] for r in rs]))
            mq = float(np.mean([r["scores"]["QSCa"][k] for r in rs]))
            bands[k] = dict(Input=mi, QSCa=mq, assessed=mi >= 0.2, ok=(mi < 0.2) or (mq >= mi - 0.05))
        h["H3_bands"] = bands; h["H3_pass"] = all(v["ok"] for v in bands.values())
        fac = [r["info"]["qscan"]["factor"] for r in rs]
        within = [abs(np.log(x / c)) <= np.log(1.25) for x in fac]
        h["H4_factors"] = fac; h["H4_within"] = int(sum(within)); h["H4_pass"] = sum(within) >= 0.8 * len(rs)
        h["H4_at_boundary"] = int(sum(r["info"]["qscan"]["at_boundary"] for r in rs))
        # H5 (descriptive): does the well-tie analogue see what the truth sees?
        m1 = [np.mean([r["scores"][m]["M1"] for r in rs]) for m in METHODS]
        t2 = [np.mean([r["scores"][m]["T2like"] for r in rs]) for m in METHODS]
        h["H5_kendall_methods_T2like_vs_M1"] = float(stats.kendalltau(m1, t2)[0])
        dm, dt = diff(rs, "QSCa", "Input"), diff(rs, "QSCa", "Input", "T2like")
        h["H5_sign_agreement_QSCa_vs_Input"] = float(np.mean(np.sign(dm) == np.sign(dt)))
        h["H5_T2like_QSCa_minus_Input"] = ci(dt)
        h["f_true_mean"] = float(np.mean([r["f_true"] for r in rs]))
        h["fc_mean"] = float(np.mean([r["info"]["wavelet"]["fc"] for r in rs]))
        H[f"c{c}"] = h
    out["primary"] = H
    out["verdict"] = dict(H1=all(H[k]["H1_pass"] for k in H), H3=all(H[k]["H3_pass"] for k in H),
                          H4=all(H[k]["H4_pass"] for k in H))
    lo = [H[k]["H2_QSCa_minus_QIQa"][1] for k in H]; hi = [H[k]["H2_QSCa_minus_QIQa"][2] for k in H]
    out["verdict"]["H2_vs_QIQa"] = ("QSCa better" if all(x > 0 for x in lo) else
                                    "QIQa better" if all(x < 0 for x in hi) else "mixed")
    lo = [H[k]["H2_QSCa_minus_IQ20a"][1] for k in H]; hi = [H[k]["H2_QSCa_minus_IQ20a"][2] for k in H]
    out["verdict"]["H2_vs_IQ20a"] = ("QSCa better" if all(x > 0 for x in lo) else
                                     "IQ20a better" if all(x < 0 for x in hi) else "mixed")
    sec = {}
    for g in ("G", "R1", "PH"):
        for c in C_TRUE:
            rs = group(rows, g, c)
            if rs:
                sec[f"{g}_c{c}"] = dict(n=len(rs), QSCa_minus_Input=ci(diff(rs, "QSCa", "Input")),
                                        QSCa_minus_QIQa=ci(diff(rs, "QSCa", "QIQa")),
                                        qscan=[r["info"]["qscan"]["factor"] for r in rs])
    out["secondary"] = sec
    bm.writejson(AN / "summary.json", out)
    print(json.dumps(out["verdict"]), flush=True)
    for k, h in H.items():
        print(k, "n", h["n"], "H1", np.round(h["H1_QSCa_minus_Input"], 3), "H2 QIQa", np.round(h["H2_QSCa_minus_QIQa"], 3),
              "H2 IQ20a", np.round(h["H2_QSCa_minus_IQ20a"], 3), "H3", h["H3_pass"], "H4", h["H4_within"], flush=True)


if __name__ == "__main__":
    main()
