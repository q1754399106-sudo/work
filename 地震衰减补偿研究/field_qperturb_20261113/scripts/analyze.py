"""Verdicts R1-R5 of the field Q-perturbation test (PLAN_ZH.md section 4) and descriptive summaries -> analysis/summary.json.
Rule as experiment 13: 90% t interval of the 8 group-mean differences.
"""
import json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; E13 = ROOT.parent / "field_baselines_volve_20261029"; T95 = {7: 1.895}
M_ALL = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); LOW, HIGH = (0.5, 0.6, 0.7), (1.4, 1.7, 2.0); METHODS = ("UQW", "VDB", "IQ70")
R = json.loads((ROOT / "results" / "tie.json").read_text(encoding="utf-8")); W, G = R["wells"], R["groups"]
t2 = lambda n, a, m: W[n]["scores"][f"{a}_{m}"]["T2"]
grp = lambda fn: np.array([np.mean([fn(n) for n in mem]) for mem in G.values()])


def stat(d):
    n = len(d); m = float(d.mean()); se = float(d.std(ddof=1) / np.sqrt(n)); k = T95[n - 1]; s = dict(mean=m, ci90=[m - k * se, m + k * se], positive=int((d > 0).sum()), values=d.tolist())
    s["verdict"] = "支持" if m > 0 and s["ci90"][0] > 0 and s["positive"] >= 6 else ("相反" if m < 0 and s["ci90"][1] < 0 else "无法判定"); return s


agg = {"low": lambda n, a: np.mean([t2(n, a, m) for m in LOW]), "all": lambda n, a: np.mean([t2(n, a, m) for m in M_ALL]),
       "worst": lambda n, a: min(t2(n, a, m) for m in M_ALL), "high": lambda n, a: np.mean([t2(n, a, m) for m in HIGH])}
TESTS = {"R1": ("low", "UQW", "VDB"), "R2": ("low", "UQW", "IQ70"), "R3": ("all", "UQW", "VDB"), "R4": ("worst", "UQW", "VDB"), "R5": ("high", "UQW", "VDB")}
out = {k: stat(grp(lambda n: agg[g](n, a) - agg[g](n, b))) for k, (g, a, b) in TESTS.items()}
out["desc_vs_IQ70"] = {g: stat(grp(lambda n: agg[g](n, "UQW") - agg[g](n, "IQ70"))) for g in ("all", "worst", "high")}
out["group_mean_T2"] = {a: {str(m): float(grp(lambda n: t2(n, a, m)).mean()) for m in M_ALL} for a in METHODS}; out["group_mean_T2"]["Input"] = float(grp(lambda n: W[n]["scores"]["Input"]["T2"]).mean())
out["diff_by_m_UQW_minus_VDB"] = {str(m): stat(grp(lambda n: t2(n, "UQW", m) - t2(n, "VDB", m))) for m in M_ALL}
E = json.loads(subprocess.run(["git", "show", "HEAD:./analysis/summary.json"], cwd=E13, capture_output=True, check=True).stdout.decode("utf-8"))["wells"]
out["check_vs_e13"] = {a: max(abs(W[n]["scores"][b]["T2"] - E[n][c]) for n in W) for a, b, c in (("Input", "Input", "Input"), ("VDB", "VDB_1.0", "B3_VDB"), ("IQ70", "IQ70_1.0", "B1_IQ70"))}
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "lines").glob("IL_*.json"))]; bk = ["10-25", "25-40", "40-55", "55-70", "70-85"]
out["n_lines"] = len(recs); out["self_corr_median"] = {a: {str(m): float(np.median([r["self_corr"][f"{a}_{m}"] for r in recs])) for m in M_ALL} for a in METHODS}
out["snr_median"] = {k: {b: float(np.median([r["bands"][k][b] for r in recs])) for b in bk} for k in ["Input"] + [f"{a}_{m}" for a in METHODS for m in (0.5, 1.0, 2.0)]}
out["band05_95_median"] = {k: np.median([r["bands"][k]["band05_95"] for r in recs], axis=0).tolist() for k in ["Input"] + [f"{a}_{m}" for a in METHODS for m in (0.5, 1.0, 2.0)]}
out["c_sr"] = [float(np.min([r["c_sr"] for r in recs])), float(np.max([r["c_sr"] for r in recs]))]
(ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print("lines", out["n_lines"], "check vs experiment 13", out["check_vs_e13"], "c_sr", out["c_sr"])
for k in TESTS:
    s = out[k]; print(k, TESTS[k], "%+.3f [%+.3f, %+.3f] %d/8 %s" % (s["mean"], *s["ci90"], s["positive"], s["verdict"]))
for g, s in out["desc_vs_IQ70"].items():
    print("UQW - IQ70", g, "%+.3f [%+.3f, %+.3f] %d/8" % (s["mean"], *s["ci90"], s["positive"]))
print("Input T2 %.3f" % out["group_mean_T2"]["Input"])
for a in METHODS:
    print("%-5s T2 by m:" % a, " ".join("%.3f" % out["group_mean_T2"][a][str(m)] for m in M_ALL), "| self corr:", " ".join("%.2f" % out["self_corr_median"][a][str(m)] for m in M_ALL))
print("UQW - VDB by m:", " ".join("%+.3f(%d)" % (out["diff_by_m_UQW_minus_VDB"][str(m)]["mean"], out["diff_by_m_UQW_minus_VDB"][str(m)]["positive"]) for m in M_ALL))
for k, v in out["snr_median"].items():
    print("%-9s SNR" % k, {b: round(x, 2) for b, x in v.items()}, "band", [round(x) for x in out["band05_95_median"][k]])
