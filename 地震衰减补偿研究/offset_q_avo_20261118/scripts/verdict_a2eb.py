"""Verdict of PLAN_GATE_A2EB_ZH.md section 3 (E1-E5 as gate A2 on the original five c, E6 shrinkage cost at c 0.85/1.2),
plus descriptive tables. Cell = (model, fp, scene, c, snr), seed mean. -> results/a2eb_verdict.json"""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; GM = ["M32", "M33", "M34", "M35"]; C5 = (0.5, 0.7, 1.0, 1.4, 2.0); CALL = (0.5, 0.7, 0.85, 1.0, 1.2, 1.4, 2.0)
R = [r for m in GM for fp in (60, 40) for r in json.loads((ROOT / "results" / ("a2eb_%s_f%d.json" % (m, fp))).read_text())]
ALN = ("events_aligned", "thin_aligned"); METH = ("oracle", "uncal", "tik_uncal", "cal", "cal_eb")


def cell(method, q="ga"):
    out = {}
    for r in R:
        if r["method"] == method:
            out.setdefault((r["model"], r["fp"], r["scene"], r["c"], r["snr"]), []).append(r[q])
    return {k: float(np.mean(v)) for k, v in out.items()}


ga = {m: cell(m) for m in METH}; wv = {m: cell(m, "wave") for m in METH}
keys = [k for k in ga["cal_eb"] if k[2] in ALN]; k5 = [k for k in keys if k[3] in C5]; out = {}
for M in ("cal_eb", "cal"):
    v = {}
    groups = sorted({(k[0], k[1], k[2], k[4]) for k in k5}); e1 = []
    for g in groups:
        ks = [(g[0], g[1], g[2], c, g[3]) for c in C5]
        e1.append(np.mean([ga[M][k] for k in ks]) / min(np.mean([ga["uncal"][k] for k in ks]), np.mean([ga["tik_uncal"][k] for k in ks])))
    v["E1_frac"] = float(np.mean(np.array(e1) <= 0.5)); v["E1_pass"] = bool(v["E1_frac"] >= 0.9)
    v["E2_median"] = float(np.median([ga[M][k] - ga["oracle"][k] for k in k5])); v["E2_pass"] = bool(v["E2_median"] <= 0.10)
    k1 = [k for k in k5 if k[3] == 1.0]
    v["E3_median"] = float(np.median([ga[M][k] - ga["uncal"][k] for k in k1])); v["E3_pass"] = bool(v["E3_median"] <= 0.05)
    v["E4_median"] = float(np.median([wv[M][k] - wv["uncal"][k] for k in k1])); v["E4_pass"] = bool(v["E4_median"] <= 0.02)
    kn = [k for k in ga[M] if k[2] == "events_nmo" and k[3] in C5 and k[3] != 1.0]
    v["E5_frac"] = float(np.mean([ga[M][k] <= ga["uncal"][k] for k in kn])); v["E5_pass"] = bool(v["E5_frac"] >= 0.8)
    out[M] = v
km = [k for k in keys if k[3] in (0.85, 1.2)]
out["cal_eb"]["E6_median"] = float(np.median([ga["cal_eb"][k] - ga["cal"][k] for k in km])); out["cal_eb"]["E6_pass"] = bool(out["cal_eb"]["E6_median"] <= 0.03)
out["all_pass"] = all(out["cal_eb"][k] for k in out["cal_eb"] if k.endswith("_pass"))
d = {}
for c in CALL:
    row = {m: dict(ga=float(np.median([ga[m][k] for k in keys if k[3] == c])), wave=float(np.median([wv[m][k] for k in keys if k[3] == c]))) for m in METH}
    rr = [r for r in R if r["method"] == "cal" and r["c"] == c and r["scene"] == "events_aligned"]
    row["rms_ln_hard"] = float(np.sqrt(np.mean([np.log((r["c_hat"] if r["c_hat"] == r["c_hat"] and r["c_hat"] > 0 else 1.0) / c) ** 2 for r in rr])))
    row["rms_ln_eb"] = float(np.sqrt(np.mean([np.log(r["c_eb"] / c) ** 2 for r in rr])))
    row["share_shrunk_to_1"] = float(np.mean([r["c_eb"] == 1.0 for r in rr]))
    d["c%.2f" % c] = row
out["descriptive"] = d
(ROOT / "results" / "a2eb_verdict.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "descriptive"}, indent=1))
for c in CALL:
    r = d["c%.2f" % c]
    print("c %.2f  " % c + "  ".join("%s %.3f" % (m, r[m]["ga"]) for m in METH) + "  ln rms hard %.3f eb %.3f  ->1 %.2f" % (r["rms_ln_hard"], r["rms_ln_eb"], r["share_shrunk_to_1"]))
