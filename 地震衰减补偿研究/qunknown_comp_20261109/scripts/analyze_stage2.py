"""Analysis of PLAN_STAGE2_ZH.md: rules P0, Q1-Q3 on the fourth batch; descriptive tables. -> analysis/stage2.json"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import M_ALL  # noqa: E402
from stage2 import IDS, SIG_B, SIG_LOGN  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; R2 = ROOT / "results" / "stage2"
rng = np.random.default_rng(20261112); NB = 2000; IDX = rng.integers(0, len(IDS), (NB, len(IDS)))
rec = {d: json.loads((R2 / f"{d}.json").read_text(encoding="utf-8"))["records"] for d in IDS}


def s(method, m, sig=None):
    out = []
    for d in IDS:
        hit = [r["S"] for r in rec[d] if r["method"] == method and r["m"] == m and (sig is None or abs(r["sig"] - sig) < 1e-12)]
        assert len(hit) == 1, (d, method, m, sig, hit); out.append(hit[0])
    return np.array(out)


def ci(x):
    b = x[IDX].mean(axis=1); lo, hi = np.quantile(b, [0.05, 0.95]); return [float(x.mean()), float(lo), float(hi)]


mean_m = lambda f: np.mean([f(m) for m in M_ALL], axis=0)                       # per dataset, uniform over m
X = mean_m(lambda m: s("UQ_lognormal", m, SIG_B))
base = {k: mean_m(lambda m, k=k: s(k, m, SIG_B if k == "C2" else None)) for k in ("IQ70", "VDB", "C1", "C2")}
out = dict(P0=dict(max_rel_diff=max(r["rel_diff"] for d in IDS for r in rec[d] if r["method"] == "UQ0")))
out["P0"]["pass"] = out["P0"]["max_rel_diff"] <= 1e-10
out["mean_S"] = dict(UQ=float(X.mean()), **{k: float(v.mean()) for k, v in base.items()})
out["Q1"] = ci(X - base["IQ70"]); out["Q2"] = ci(X - base["VDB"]); out["Q3_C1"] = ci(X - base["C1"]); out["Q3_C2"] = ci(X - base["C2"])
out["pass"] = dict(Q1=out["Q1"][0] >= 0.03 and out["Q1"][1] > 0, Q2=out["Q2"][1] > 0,
                   Q3=out["Q3_C1"][1] > 0 and out["Q3_C2"][1] > 0)
out["confirmed"] = out["P0"]["pass"] and all(out["pass"].values())
# descriptive
out["cost_at_m1"] = ci(s("UQ_lognormal", 1.0, SIG_B) - s("IQ70", 1.0))
units = np.concatenate([s("UQ_lognormal", m, SIG_B) - s("IQ70", m) for m in M_ALL]); out["win_rate_vs_IQ70"] = float((units > 0).mean())
out["by_structure"] = {k: {n: float((X - base[n])[[i for i, d in enumerate(IDS) if d[1] == k]].mean()) for n in ("IQ70", "VDB", "C2")}
                       for k in "JKL"}
out["declared"] = {f"{q:g}": float(mean_m(lambda m: s("UQ_lognormal", m, q)).mean()) for q in SIG_LOGN}
out["declared"]["gauss_0.447"] = float(mean_m(lambda m: s("UQ_gauss", m, SIG_B)).mean())
cur = {f"UQ_lognormal_{q:g}": ("UQ_lognormal", q) for q in SIG_LOGN}
cur.update(UQ_gauss_0447=("UQ_gauss", SIG_B), IQ70=("IQ70", None), VDB=("VDB", None), C1=("C1", None), C2=("C2", SIG_B))
out["curves"] = {n: {f"{m:g}": float(s(k, m, q).mean()) for m in M_ALL} for n, (k, q) in cur.items()}
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "stage2.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("P0", out["P0"]); print("mean S", {k: round(v, 4) for k, v in out["mean_S"].items()})
for k in ("Q1", "Q2", "Q3_C1", "Q3_C2", "cost_at_m1"):
    print(k, [round(x, 4) for x in out[k]])
print("pass", out["pass"], "confirmed", out["confirmed"], "win", round(out["win_rate_vs_IQ70"], 3))
print("by structure", out["by_structure"]); print("declared", out["declared"])
for n, c in out["curves"].items():
    print("%-18s" % n, " ".join("%.3f" % v for v in c.values()))
