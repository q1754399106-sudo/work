"""Verdicts F1-F5 of the UQW-S field test (PLAN_ZH.md) and descriptive summaries -> analysis/summary.json.
Rule as experiment 13 and field_qperturb_20261113: 90% t interval of the 8 group-mean differences; "支持" if mean > 0, lower
bound > 0 and >= 6/8 positive; "相反" if mean < 0 and upper bound < 0; otherwise "无法判定".
Consistency check: UQW, VDB and IQ70 ties at every m must equal field_qperturb_20261113/results/tie.json.
usage: analyze.py [TIE_JSON]   (argument: another tie file, for the code check on fake data; then nothing is written)"""
import sys, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; QP = ROOT.parent / "field_qperturb_20261113"; T95 = {7: 1.895}
M_ALL = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); LOW, NEAR, HIGH = (0.5, 0.6, 0.7), (0.85, 1.0, 1.2), (1.4, 1.7, 2.0)
METHODS = ("UQWS", "UQW", "VDB", "VDBin", "HWin", "IQ70")
SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "tie.json"
R = json.loads(SRC.read_text(encoding="utf-8")); W, G = R["wells"], R["groups"]
t2 = lambda n, a, m: W[n]["scores"][f"{a}_{m}"]["T2"]
grp = lambda fn: np.array([np.mean([fn(n) for n in mem]) for mem in G.values()])


def stat(d):
    n = len(d); m = float(d.mean()); se = float(d.std(ddof=1) / np.sqrt(n)); k = T95[n - 1]; s = dict(mean=m, ci90=[m - k * se, m + k * se], positive=int((d > 0).sum()), values=d.tolist())
    s["verdict"] = "支持" if m > 0 and s["ci90"][0] > 0 and s["positive"] >= 6 else ("相反" if m < 0 and s["ci90"][1] < 0 else "无法判定"); return s


agg = {"low": lambda n, a: np.mean([t2(n, a, m) for m in LOW]), "near": lambda n, a: np.mean([t2(n, a, m) for m in NEAR]),
       "high": lambda n, a: np.mean([t2(n, a, m) for m in HIGH]), "all": lambda n, a: np.mean([t2(n, a, m) for m in M_ALL]),
       "worst": lambda n, a: min(t2(n, a, m) for m in M_ALL)}
TESTS = {"F1": ("high", "UQWS", "UQW"), "F2": ("all", "UQWS", "UQW"), "F3": ("high", "UQWS", "VDB"), "F4": ("all", "UQWS", "VDB"), "F5": ("worst", "UQWS", "VDB")}
out = {k: dict(test=f"{g}: {a} - {b}", **stat(grp(lambda n: agg[g](n, a) - agg[g](n, b)))) for k, (g, a, b) in TESTS.items()}
out["describe"] = {f"{g} {a}-{b}": stat(grp(lambda n: agg[g](n, a) - agg[g](n, b)))
                   for g in ("low", "near", "high", "all", "worst") for a, b in (("UQWS", "UQW"), ("UQWS", "VDB"), ("UQWS", "VDBin"), ("UQW", "VDB"), ("VDBin", "VDB"), ("UQWS", "IQ70"))}
out["group_mean_T2"] = {a: {str(m): float(grp(lambda n: t2(n, a, m)).mean()) for m in M_ALL} for a in METHODS}
out["group_mean_T2"]["Input"] = float(grp(lambda n: W[n]["scores"]["Input"]["T2"]).mean())
if len(sys.argv) == 1:
    Q = json.loads((QP / "results" / "tie.json").read_text(encoding="utf-8"))["wells"]
    out["check_vs_qperturb"] = {a: max(abs(t2(n, a, m) - Q[n]["scores"][f"{a}_{m}"]["T2"]) for n in W for m in M_ALL) for a in ("UQW", "VDB", "IQ70")}
    recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "lines").glob("IL_*.json"))]
    out["n_lines"] = len(recs); out["fs_r"] = dict(min=float(min(r["fs_r"] for r in recs)), median=float(np.median([r["fs_r"] for r in recs])), n_le_72=int(sum(r["fs_r"] <= 72 for r in recs)))
    out["self_corr_median"] = {a: {str(m): float(np.median([r["self_corr"][f"{a}_{m}"] for r in recs])) for m in M_ALL} for a in METHODS}
    out["band05_95_median"] = {k: np.median([r["bands"][k]["band05_95"] for r in recs], axis=0).tolist() for k in ["Input"] + [f"{a}_{m}" for a in METHODS for m in (0.5, 1.0, 2.0)]}
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print("lines", out["n_lines"], "check vs qperturb", out["check_vs_qperturb"], "fs_r", out["fs_r"])
for k in TESTS:
    s = out[k]; print(k, "%-22s %+.3f [%+.3f, %+.3f] %d/8 %s" % (s["test"], s["mean"], *s["ci90"], s["positive"], s["verdict"]))
print("Input T2 %.3f" % out["group_mean_T2"]["Input"])
for a in METHODS:
    print("%-6s T2 by m:" % a, " ".join("%.3f" % out["group_mean_T2"][a][str(m)] for m in M_ALL))
for k, s in out["describe"].items():
    print("  %-18s %+.3f [%+.3f, %+.3f] %d/8 %s" % (k, s["mean"], *s["ci90"], s["positive"], s["verdict"]))
