"""Analysis of PLAN_STAGE1_ZH.md: rules P0-P4 per scenario (A, B) and law (gauss, lognormal); descriptive tables.
-> analysis/stage1.json"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import DEV, M_ALL, M_SC  # noqa: E402
from stage1 import SIG, SIG_A, SIG_B  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
R0, R1 = ROOT / "results" / "stage0", ROOT / "results" / "stage1"
SCEN = dict(A=dict(ms=(0.7, 0.85, 1.0, 1.2, 1.4), sig=SIG_A), B=dict(ms=M_ALL, sig=SIG_B))
rng = np.random.default_rng(20261111); NB = 2000; IDX = rng.integers(0, len(DEV), (NB, len(DEV)))
rec = {d: json.loads((R1 / f"{d}.json").read_text(encoding="utf-8"))["records"] for d in DEV}


def s1(method, m, sig=None):
    """S per dataset from stage 1 (sig None: the record without sig, or the only one)."""
    out = []
    for d in DEV:
        hit = [r["S"] for r in rec[d] if r["method"] == method and r["m"] == m and (sig is None or abs(r["sig"] - sig) < 1e-12)]
        assert len(hit) == 1, (d, method, m, sig, hit); out.append(hit[0])
    return np.array(out)


def s0(method, m):
    return np.array([json.loads((R0 / f"{d}__{method}_m{m:g}.json").read_text())["score"]["S"] for d in DEV])


def ci(x):
    b = x[IDX].mean(axis=1); lo, hi = np.quantile(b, [0.05, 0.95]); return [float(x.mean()), float(lo), float(hi)]


def smean(f, ms):
    return np.mean([f(m) for m in ms], axis=0)                                   # per dataset, uniform over the grid


out = dict(P0=dict(max_rel_diff=max(r["rel_diff"] for d in DEV for r in rec[d] if r["method"] == "UQ0"),
                   max_abs_dS=float(max(np.abs(s1("UQ0", m) - s0("IQ70", m)).max() for m in M_ALL))))
out["P0"]["pass"] = out["P0"]["max_rel_diff"] <= 1e-10 and out["P0"]["max_abs_dS"] <= 1e-9
out["rules"], out["curves"] = {}, {}
for sc, cfg in SCEN.items():
    ms, sig = cfg["ms"], cfg["sig"]
    for law in ("gauss", "lognormal"):
        X = smean(lambda m: s1(f"UQ_{law}", m, sig), ms)
        r = dict(mean_S=float(X.mean()),
                 P1=ci(X - smean(lambda m: s0("IQ70", m), ms)),
                 P2=float((s1(f"UQ_{law}", 1.0, sig) - s0("IQ70", 1.0)).mean()),
                 P3=ci(X - smean(lambda m: s0("VDB", m), ms)),
                 P4_C1=ci(X - smean(lambda m: s1("C1", m), ms)),
                 P4_C2=ci(X - smean(lambda m: s1("C2", m, sig), ms)))
        r["P1_pass"] = r["P1"][0] >= 0.03 and r["P1"][1] > 0
        r["P2_pass"] = r["P2"] >= -0.01
        r["P3_pass"] = r["P3"][1] > 0
        r["P4_pass"] = r["P4_C1"][1] > 0 and r["P4_C2"][1] > 0
        r["all_pass"] = r["P1_pass"] and r["P3_pass"] and r["P4_pass"] and (r["P2_pass"] or sc == "B")   # P2 binds in A only
        # descriptive: unit win rate against IQ70, SC-a on the SC grid, mis-declared sigma
        units = np.concatenate([s1(f"UQ_{law}", m, sig) - s0("IQ70", m) for m in ms])
        r["win_rate_vs_IQ70"] = float((units > 0).mean())
        msc = [m for m in ms if m in M_SC]
        r["vs_SCa_on_SC_grid"] = ci(smean(lambda m: s1(f"UQ_{law}", m, sig), msc) - smean(lambda m: s0("SCa", m), msc))
        r["declared"] = {f"{q:g}": float(smean(lambda m: s1(f"UQ_{law}", m, q), ms).mean()) for q in SIG}
        out["rules"][f"{sc}_{law}"] = r
    out["rules"][f"{sc}_baselines"] = {k: float(smean(f, ms).mean()) for k, f in
                                       (("IQ70", lambda m: s0("IQ70", m)), ("VDB", lambda m: s0("VDB", m)),
                                        ("C1", lambda m: s1("C1", m)), ("C2", lambda m: s1("C2", m, sig)),
                                        ("gauss_abs", lambda m: s1("UQ_gauss_abs", m, sig)))}
for law in ("gauss", "lognormal"):
    out["curves"][law] = {f"{q:g}": {f"{m:g}": float(s1(f"UQ_{law}", m, q).mean()) for m in M_ALL} for q in SIG}
out["curves"]["C1"] = {f"{m:g}": float(s1("C1", m).mean()) for m in M_ALL}
for q in (SIG_A, SIG_B):
    out["curves"][f"C2_{q:g}"] = {f"{m:g}": float(s1("C2", m, q).mean()) for m in M_ALL}
    out["curves"][f"gauss_abs_{q:g}"] = {f"{m:g}": float(s1("UQ_gauss_abs", m, q).mean()) for m in M_ALL}
sr = {k: np.array([[r["S"] for r in rec[d] if r["method"] == k][0] for d in DEV])
      for k in ("SR_IQ70", "SR_VDB", "SR_UQ_gauss", "SR_UQ_lognormal")}
dev = {d: json.loads((BASE / "mu_tune_20261031" / "results" / "dev" / f"{d}.json").read_text()) for d in DEV}
sr["SR_SCa"] = np.array([dev[d]["scores"]["SR_mu0.25"]["S"] for d in DEV])
out["SR_route"] = {k: float(v.mean()) for k, v in sr.items()}
out["SR_route"]["UQ_gauss_minus_IQ70"] = ci(sr["SR_UQ_gauss"] - sr["SR_IQ70"])
out["SR_route"]["UQ_lognormal_minus_IQ70"] = ci(sr["SR_UQ_lognormal"] - sr["SR_IQ70"])
out["SR_route"]["IQ70_vs_exp17_max_abs"] = float(np.abs(sr["SR_IQ70"] - np.array([dev[d]["scores"]["IQ70"]["S"] for d in DEV])).max())
anyall = [k for k, v in out["rules"].items() if isinstance(v, dict) and v.get("all_pass")]
sets = {k: v for k, v in out["rules"].items() if "P1_pass" in v}
cons = [k for k, v in sets.items() if v["P1_pass"] and v["P3_pass"] and not v["P4_pass"]]
trade = [k for k, v in sets.items() if v["P1_pass"] and v["P3_pass"] and v["P4_pass"] and not v["all_pass"]]   # P2 fails in A
out["decision"] = ("invalid" if not out["P0"]["pass"] else "continue" if anyall else "tradeoff" if trade
                   else "conservative_only" if cons else "stop")
out["passing"], out["tradeoff_sets"], out["conservative_sets"] = anyall, trade, cons
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "stage1.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("P0", out["P0"])
for k, v in out["rules"].items():
    print(k, {a: (round(b, 4) if isinstance(b, float) else [round(x, 4) for x in b] if isinstance(b, list) else b)
              for a, b in v.items() if a != "declared"})
print("SR", out["SR_route"]); print("decision", out["decision"], anyall)
