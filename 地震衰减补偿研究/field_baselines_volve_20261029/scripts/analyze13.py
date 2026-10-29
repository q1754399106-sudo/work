"""Verdicts of experiment 13 (PLAN_ZH.md section 4) and the descriptive summaries -> analysis/summary.json."""
import sys, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
T95 = {5: 2.015, 6: 1.943, 7: 1.895}             # t(0.95, n-1), as experiment 9
COMPARE = {"F1": ("B1_IQ70", "PROP"), "F2": ("B2_SCa", "PROP"), "F3": ("QSCa_tuned", "QSCa"), "F4": ("B1_IQ70", "IQ30")}
DESC = [("B2_SCa", "B1_IQ70"), ("B3_VDB", "B1_IQ70"), ("B2_SCa", "QSCa_tuned"), ("PROP", "QSCa"), ("B1_IQ70", "QSCa"),
        ("B2_SCa", "QSCa")]


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


out = dict(groups=G, wells={n: {m: v["T2"] for m, v in r["scores"].items()} for n, r in W.items()},
           shifts={n: {m: v["shift_ms"] for m, v in r["scores"].items()} for n, r in W.items()})
out["group_mean_T2"] = {m: float(np.mean([np.mean([W[n]["scores"][m]["T2"] for n in mem]) for mem in G.values()]))
                        for m in names}
for key, (a, b) in COMPARE.items():
    dg = by_group(a, b); s = stat(list(dg.values())); s["by_group"] = dg; s["verdict"] = verdict(s); out[key] = s
out["vs_input"] = {m: stat(list(by_group(m, "Input").values())) for m in names if m != "Input"}
out["descriptive"] = {f"{a} - {b}": stat(list(by_group(a, b).values())) for a, b in DESC}
E9 = json.loads((ROOT / "cache" / "inputs" / "multiwell_tie_volve_20261028" / "results" / "tie.json").read_text(encoding="utf-8"))
out["check_input_vs_e9"] = max(abs(W[n]["scores"]["Input"]["T2"] - E9["wells"][n]["scores"]["Input"]["T2"]) for n in W)
recs = {}
for p in sorted((ROOT / "results" / "lines").glob("IL_*_*.json")):
    r = json.loads(p.read_text(encoding="utf-8")); recs.setdefault(r["inline"], {}).update(r["bands"])
bkeys = ["10-25", "25-40", "40-55", "55-70", "70-85"]
out["bands_all"] = {m: {b: float(np.median([recs[ln][m][b] for ln in recs])) for b in bkeys} for m in names}
out["bands_10141"] = recs[10141]
out["qscan"] = {p.stem: json.loads(p.read_text(encoding="utf-8"))["c"] for p in sorted((ROOT / "results" / "qscan").glob("*.json"))}
(ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print({k: out[k]["verdict"] for k in COMPARE}, "input check", out["check_input_vs_e9"])
for k in COMPARE:
    s = out[k]; print(k, COMPARE[k], f"{s['mean']:+.3f} [{s['ci90'][0]:+.3f}, {s['ci90'][1]:+.3f}] {s['positive']}/8")
for m in sorted(names, key=lambda m: -out["group_mean_T2"][m]):
    print(f"{m:12s} T2 {out['group_mean_T2'][m]:.3f}  SNR", {b: round(v, 2) for b, v in out["bands_all"][m].items()})
for k, s in out["descriptive"].items():
    print(k, f"{s['mean']:+.3f} [{s['ci90'][0]:+.3f}, {s['ci90'][1]:+.3f}] {s['positive']}/8")
print("qscan", out["qscan"])
