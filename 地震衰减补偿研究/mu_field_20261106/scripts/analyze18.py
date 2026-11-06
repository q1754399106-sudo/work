"""Verdicts of experiment 18 (PLAN_ZH.md), and the descriptive summaries -> analysis/summary.json."""
import json
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RES = ROOT.parent; REPO = RES.parent; COMMIT = "27a065a"; E13 = "field_baselines_volve_20261029"
T95 = {5: 2.015, 6: 1.943, 7: 1.895}             # t(0.95, n-1), as experiments 9 and 13
COMPARE = {"M1": ("SR_mu025", "B2_SCa"), "M2": ("QT_mu025", "QSCa_tuned"), "M3": ("QT_mu025", "QSCa")}
BANDS = ["10-25", "25-40", "40-55", "55-70", "70-85"]


def git_json(rel):
    return json.loads(subprocess.check_output(["git", "show", f"{COMMIT}:{RES.name}/{rel}"], cwd=REPO).decode("utf-8"))


def stat(d):
    d = np.asarray(d, float); n = len(d); m = float(d.mean()); se = float(d.std(ddof=1) / np.sqrt(n)); k = T95[n - 1]
    return dict(n=n, mean=m, ci90=[m - k * se, m + k * se], positive=int((d > 0).sum()), values=d.tolist())


def verdict(s):
    if s["mean"] > 0 and s["ci90"][0] > 0 and s["positive"] >= 6:
        return "支持"
    if s["mean"] < 0 and s["ci90"][1] < 0:
        return "相反"
    return "无法判定"


R = json.loads((ROOT / "results" / "tie.json").read_text(encoding="utf-8")); W = R["wells"]; G = R["groups"]
names = list(next(iter(W.values()))["scores"])


def by_group(a, b):
    return {g: float(np.mean([W[n]["scores"][a]["T2"] - W[n]["scores"][b]["T2"] for n in mem])) for g, mem in G.items()}


out = dict(groups=G, check_vs_e13=R["check_vs_e13"], wells={n: {m: v["T2"] for m, v in r["scores"].items()} for n, r in W.items()},
           shifts={n: {m: v["shift_ms"] for m, v in r["scores"].items()} for n, r in W.items()})
out["group_mean_T2"] = {m: float(np.mean([np.mean([W[n]["scores"][m]["T2"] for n in mem]) for mem in G.values()])) for m in names}
for key, (a, b) in COMPARE.items():
    dg = by_group(a, b); s = stat(list(dg.values())); s["by_group"] = dg; s["verdict"] = verdict(s); out[key] = s
out["vs_input"] = {m: stat(list(by_group(m, "Input").values())) for m in names if m != "Input"}

v = [out[k]["verdict"] for k in ("M1", "M2")]
out["default_rule"] = "keep 10, conflict recorded" if "相反" in v else "change the default to 0.25"

new, old = {}, {}
for p in sorted((ROOT / "results" / "lines").glob("IL_*.json")):
    r = json.loads(p.read_text(encoding="utf-8")); new[r["inline"]] = r["bands"]
    a = git_json(f"{E13}/results/lines/IL_{r['inline']}_A.json"); b = git_json(f"{E13}/results/lines/IL_{r['inline']}_B.json")
    old[r["inline"]] = {**a["bands"], **b["bands"]}
    assert abs(a["c_sr"] - r["c_sr"]) < 1e-9 and abs(b["c_qscan"] - r["c_qscan"]) < 1e-9, r["inline"]
med = lambda src, m, b: float(np.median([src[ln][m][b] for ln in src]))
out["bands"] = {**{m: {b: med(new, m, b) for b in BANDS} for m in ("SR_mu025", "QT_mu025")},
                **{m: {b: med(old, m, b) for b in BANDS} for m in ("Input", "B2_SCa", "QSCa_tuned", "QSCa")}}
out["band05_95"] = {**{m: np.median([new[ln][m]["band05_95"] for ln in new], axis=0).tolist() for m in ("SR_mu025", "QT_mu025")},
                    **{m: np.median([old[ln][m]["band05_95"] for ln in old], axis=0).tolist() for m in ("B2_SCa", "QSCa_tuned", "QSCa")}}
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print("check vs experiment 13:", R["check_vs_e13"], "default rule:", out["default_rule"])
for k, (a, b) in COMPARE.items():
    s = out[k]; print(k, a, "-", b, f"{s['mean']:+.3f} [{s['ci90'][0]:+.3f}, {s['ci90'][1]:+.3f}] {s['positive']}/8 {s['verdict']}", {g: round(x, 3) for g, x in s["by_group"].items()})
for m in sorted(names, key=lambda m: -out["group_mean_T2"][m]):
    extra = out["vs_input"].get(m, {}); print(f"{m:11s} T2 {out['group_mean_T2'][m]:.3f}", f"vs Input {extra['mean']:+.3f} {extra['positive']}/8" if extra else "", {b: round(x, 2) for b, x in out["bands"][m].items()})
print("band 5-95%:", {m: [round(x, 1) for x in v] for m, v in out["band05_95"].items()})
