"""Verdict of PLAN_GATE_A2_ZH.md section 4, and descriptive tables. Cell = (model, fp, scene, c, snr), seed mean."""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; GM = ["M24", "M25", "M26", "M27"]; CS = (0.5, 0.7, 1.0, 1.4, 2.0)
R = [r for m in GM for fp in (60, 40) for r in json.loads((ROOT / "results" / ("a2_%s_f%d.json" % (m, fp))).read_text())]
ALN = ("events_aligned", "thin_aligned")


def cell(method, q="ga"):
    out = {}
    for r in R:
        if r["method"] == method:
            out.setdefault((r["model"], r["fp"], r["scene"], r["c"], r["snr"]), []).append(r[q])
    return {k: float(np.mean(v)) for k, v in out.items()}


ga = {m: cell(m) for m in ("oracle", "uncal", "tik_uncal", "cal")}; wv = {m: cell(m, "wave") for m in ("oracle", "uncal", "tik_uncal", "cal")}
keys = [k for k in ga["cal"] if k[2] in ALN]; out = {}
# E1: per group (model, fp, scene, snr), average over the 5 c: cal <= 0.5 x min(uncal, tik_uncal) in >= 90 % of groups
groups = sorted({(k[0], k[1], k[2], k[4]) for k in keys}); e1 = []
for g in groups:
    ks = [(g[0], g[1], g[2], c, g[3]) for c in CS]; cal = np.mean([ga["cal"][k] for k in ks])
    base = min(np.mean([ga["uncal"][k] for k in ks]), np.mean([ga["tik_uncal"][k] for k in ks])); e1.append(cal / base)
out["E1_frac"] = float(np.mean(np.array(e1) <= 0.5)); out["E1_median_ratio"] = float(np.median(e1)); out["E1_pass"] = bool(out["E1_frac"] >= 0.9)
d2 = [ga["cal"][k] - ga["oracle"][k] for k in keys]
out["E2_median_cal_minus_oracle"] = float(np.median(d2)); out["E2_pass"] = bool(out["E2_median_cal_minus_oracle"] <= 0.10)
k1 = [k for k in keys if k[3] == 1.0]
out["E3_median_cost_c1"] = float(np.median([ga["cal"][k] - ga["uncal"][k] for k in k1])); out["E3_pass"] = bool(out["E3_median_cost_c1"] <= 0.05)
out["E4_median_wave_cost_c1"] = float(np.median([wv["cal"][k] - wv["uncal"][k] for k in k1])); out["E4_pass"] = bool(out["E4_median_wave_cost_c1"] <= 0.02)
kn = [k for k in ga["cal"] if k[2] == "events_nmo" and k[3] != 1.0]
out["E5_frac_nmo_cal_le_uncal"] = float(np.mean([ga["cal"][k] <= ga["uncal"][k] for k in kn])); out["E5_pass"] = bool(out["E5_frac_nmo_cal_le_uncal"] >= 0.8)
out["all_pass"] = all(out[k] for k in out if k.endswith("_pass"))
# descriptive
d = {}
for c in CS:
    row = {}
    for m in ("oracle", "uncal", "tik_uncal", "cal"):
        row[m] = dict(ga=float(np.median([ga[m][k] for k in keys if k[3] == c])), wave=float(np.median([wv[m][k] for k in keys if k[3] == c])))
    ch = [np.log(r["c_hat"] / r["c"]) for r in R if r["method"] == "cal" and r["c"] == c and r["scene"] == "events_aligned" and r["c_hat"] == r["c_hat"] and r["c_hat"] > 0]
    row["ln_chat_err"] = dict(mean=float(np.mean(ch)), sd=float(np.std(ch)), rms=float(np.sqrt(np.mean(np.square(ch)))))
    row["chat_fail"] = int(sum(1 for r in R if r["method"] == "cal" and r["c"] == c and r["scene"] == "events_aligned" and not (r["c_hat"] == r["c_hat"] and r["c_hat"] > 0)))
    d["c%.1f" % c] = row
for grp, idx in (("fp", 1), ("scene", 2), ("snr", 4), ("model", 0)):
    g = {}
    for k in keys:
        g.setdefault(str(k[idx]), {m: [] for m in ga})
        for m in ga:
            g[str(k[idx])][m].append(ga[m][k])
    d["by_" + grp] = {a: {m: float(np.median(v)) for m, v in b.items()} for a, b in g.items()}
d["nmo"] = {m: float(np.median([ga[m][k] for k in ga["cal"] if k[2] == "events_nmo"])) for m in ga}
out["descriptive"] = d
(ROOT / "results" / "a2_verdict.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "descriptive"}, indent=1))
