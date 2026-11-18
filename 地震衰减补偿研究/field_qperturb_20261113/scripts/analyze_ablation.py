"""Verdict of the post hoc ablation (PLAN_ABLATION_ZH.md section 4) -> analysis/ablation.json. Group means over the 8 well
groups and 90% t intervals, as analyze.py.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; T95 = 1.895
M_ALL = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); SIDES = {"mid": (0.85, 1.0, 1.2), "low": (0.5, 0.6, 0.7), "high": (1.4, 1.7, 2.0), "all": M_ALL}
A = json.loads((ROOT / "results" / "ablation_tie.json").read_text(encoding="utf-8")); W, G = A["wells"], A["groups"]
OLD = json.loads((ROOT / "results" / "tie.json").read_text(encoding="utf-8"))["wells"]
grp = lambda fn: np.array([np.mean([fn(n) for n in mem]) for mem in G.values()])


def stat(d):
    m = float(d.mean()); se = float(d.std(ddof=1) / np.sqrt(len(d))); return dict(mean=m, ci90=[m - T95 * se, m + T95 * se], positive=int((d > 0).sum()), values=d.tolist())


def gap(ref, a, b, side):
    return stat(grp(lambda n: np.mean([W[n][ref][f"{a}_{m}"] - W[n][ref][f"{b}_{m}"] for m in SIDES[side]])))


out = {"check_vs_tie_json": max(abs(W[n]["ricker"][f"{a}_{m}"] - OLD[n]["scores"][f"{a}_{m}"]["T2"]) for n in W for a in ("UQW", "VDB", "IQ70") for m in M_ALL)}
out["check_vdb_rebuild"] = max(abs(W[n]["ricker"][f"VDB_re_{m}"] - W[n]["ricker"][f"VDB_{m}"]) for n in W for m in M_ALL)
VARS = ("UQW", "UQW_TI", "UQW_TId", "UQW_TV", "PH", "PH_TI", "IQ70")
for ref in ("ricker", "wout"):
    out[ref] = {s: {v: gap(ref, v, "VDB", s) for v in VARS} for s in SIDES}
    out[ref]["group_mean_T2"] = {v: {str(m): float(grp(lambda n: W[n][ref][f"{v}_{m}"]).mean()) for m in M_ALL} for v in VARS + ("VDB",)}

# primary verdict: Ricker reference, Q near correct (m 0.85, 1.0, 1.2)
base = out["ricker"]["mid"]["UQW"]; ti = out["ricker"]["mid"]["UQW_TI"]; tv = out["ricker"]["mid"]["UQW_TV"]
red = lambda s: 1 - s["mean"] / base["mean"]
v = dict(base=base["mean"], base_ci90=base["ci90"], r_TI=red(ti), r_TV=red(tv), r_TId=red(out["ricker"]["mid"]["UQW_TId"]), TI_ci90=ti["ci90"])
if not base["ci90"][1] < 0:
    v["verdict"] = "无法判定（基准差距的 90% 区间含 0）"
elif v["r_TI"] >= 0.6 and ti["ci90"][1] >= 0:
    v["verdict"] = "静态为主（解释 B）"
elif v["r_TI"] < 0.3 and v["r_TV"] >= 0.6:
    v["verdict"] = "随深度变化为主（解释 A）"
else:
    v["verdict"] = "混合 / 无法判定"
out["verdict"] = v
desc = ROOT / "results" / "ablation_desc.json"
if desc.exists():
    out["desc"] = json.loads(desc.read_text(encoding="utf-8"))
(ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "ablation.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

print("checks: tie.json max |diff| %.2e, VDB rebuild tie max |diff| %.2e" % (out["check_vs_tie_json"], out["check_vdb_rebuild"]))
print("VERDICT:", v["verdict"], "| base %+.3f [%+.3f, %+.3f]  r_TI %.2f  r_TId %.2f  r_TV %.2f" % (v["base"], *v["base_ci90"], v["r_TI"], v["r_TId"], v["r_TV"]))
for ref in ("ricker", "wout"):
    print("reference", ref)
    for s in SIDES:
        print("  %-4s" % s, "  ".join("%s %+.3f(%d)" % (k, out[ref][s][k]["mean"], out[ref][s][k]["positive"]) for k in VARS))
