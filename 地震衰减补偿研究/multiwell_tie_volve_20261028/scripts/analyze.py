"""H1 of experiment 9 (PLAN_ZH.md section 4) and the descriptive comparisons -> analysis/summary.json."""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq  # noqa: E402
from wq import np  # noqa: E402

T95 = {5: 2.015, 6: 1.943, 7: 1.895}


def stat(d):
    d = np.asarray(d, float); n = len(d); m = float(d.mean()); se = float(d.std(ddof=1) / np.sqrt(n)); k = T95[n - 1]
    return dict(n=n, mean=m, ci90=[m - k * se, m + k * se], positive=int((d > 0).sum()), values=d.tolist())


def verdict(s):
    if s["mean"] > 0 and s["ci90"][0] > 0 and s["positive"] >= 6:
        return "成立"
    if s["mean"] < 0 and s["ci90"][1] < 0:
        return "相反"
    return "无法判定"


R = json.loads((wq.ROOT / "results" / "tie.json").read_text(encoding="utf-8")); W = R["wells"]; G = R["groups"]
out = dict(groups=G, wells={n: {k: v["T2"] for k, v in r["scores"].items()} for n, r in W.items()})
for m in ("IQ30_SR", "IQ30_0.77", "IQ30_1.0"):
    dg = {g: float(np.mean([W[n]["scores"][m]["T2"] - W[n]["scores"]["Input"]["T2"] for n in mem])) for g, mem in G.items()}
    s = stat(list(dg.values())); s["by_group"] = dg
    unseen = [v for g, v in dg.items() if not any(n in R["seen"] for n in G[g])]
    s["unseen"] = stat(unseen)
    out[f"{m}_minus_Input"] = s
out["H1"] = verdict(out["IQ30_SR_minus_Input"])
for m in ("IQ30_0.77", "IQ30_1.0"):
    dg = [float(np.mean([W[n]["scores"][m]["T2"] - W[n]["scores"]["IQ30_SR"]["T2"] for n in mem])) for mem in G.values()]
    out[f"{m}_minus_IQ30_SR"] = stat(dg)
(wq.ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print(json.dumps({k: v for k, v in out.items() if k not in ("wells",)}, ensure_ascii=False, indent=1))
