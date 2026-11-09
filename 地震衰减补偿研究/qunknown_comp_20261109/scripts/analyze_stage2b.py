"""Analysis of PLAN_STAGE2B_ZH.md: rules P0, G1-G3 on the FD models (14 dB); descriptive tables. -> analysis/stage2b.json"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import M_ALL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; RD = ROOT / "results" / "stage2b"
MODELS = [f"H{i}" for i in range(1, 7)]; R = (1, 2, 3); SIG = 0.447
rng = np.random.default_rng(20261113); NB = 2000; IDX = rng.integers(0, 6, (NB, 6))
rec = {(h, lv, r): json.loads((RD / f"{h}_{lv}_r{r}.json").read_text(encoding="utf-8"))["records"]
       for h in MODELS for lv in ("14dB", "08dB") for r in R}


def val(h, lv, method, m, sig=None):
    """Mean over the 3 realisations for one model, method, m."""
    v = []
    for r in R:
        hit = [x["S"] for x in rec[(h, lv, r)] if x["method"] == method and x["m"] == m
               and (sig is None or abs(x["sig"] - sig) < 1e-12)]
        assert len(hit) == 1, (h, lv, r, method, m, sig, hit); v.append(hit[0])
    return float(np.mean(v))


def per_model(lv, method, sig=None, ms=M_ALL):
    return np.array([np.mean([val(h, lv, method, m, sig) for m in ms]) for h in MODELS])


def ci(x):
    b = x[IDX].mean(axis=1); lo, hi = np.quantile(b, [0.05, 0.95]); return [float(x.mean()), float(lo), float(hi)]


out = dict(P0=dict(max_rel_diff=max(x["rel_diff"] for k in rec for x in rec[k] if x["method"] == "UQ0")))
out["P0"]["pass"] = out["P0"]["max_rel_diff"] <= 1e-10
for lv in ("14dB", "08dB"):
    U = per_model(lv, "UQ", SIG)
    B = {"IQ70": per_model(lv, "IQ70"), "VDB": per_model(lv, "VDB"), "C1": per_model(lv, "C1"), "C2": per_model(lv, "C2", SIG)}
    d = {k: U - v for k, v in B.items()}
    blk = dict(mean_S=dict(UQ=float(U.mean()), **{k: float(v.mean()) for k, v in B.items()}),
               diff={k: dict(ci=ci(v), per_model=dict(zip(MODELS, map(float, v))), n_pos=int((v > 0).sum())) for k, v in d.items()},
               cost_at_m1=ci(np.array([val(h, lv, "UQ", 1.0, SIG) - val(h, lv, "IQ70", 1.0) for h in MODELS])),
               declared={f"{q:g}": float(per_model(lv, "UQ", q).mean()) for q in (0.245, 0.447, 0.894)},
               anchor09=float(per_model(lv, "UQ_anchor09", SIG).mean()),
               curves={n: {f"{m:g}": float(np.mean([val(h, lv, k, m, q) for h in MODELS])) for m in M_ALL}
                       for n, (k, q) in dict(UQ_0447=("UQ", SIG), UQ_0245=("UQ", 0.245), UQ_0894=("UQ", 0.894),
                                             UQ_anchor09=("UQ_anchor09", SIG), IQ70=("IQ70", None), VDB=("VDB", None),
                                             C1=("C1", None), C2=("C2", SIG)).items()})
    if lv == "14dB":
        blk["rules"] = dict(G1=d["IQ70"].mean() >= 0.03 and (d["IQ70"] > 0).sum() >= 5,
                            G2=d["VDB"].mean() >= -0.01,
                            G3=d["C1"].mean() > 0 and (d["C1"] > 0).sum() >= 5 and d["C2"].mean() > 0 and (d["C2"] > 0).sum() >= 5)
        blk["rules"] = {k: bool(v) for k, v in blk["rules"].items()}
        blk["core_generator_independent"] = out["P0"]["pass"] and blk["rules"]["G1"] and blk["rules"]["G3"]
    out[lv] = blk
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "stage2b.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("P0", out["P0"])
for lv in ("14dB", "08dB"):
    b = out[lv]; print("==", lv, "mean S", {k: round(v, 4) for k, v in b["mean_S"].items()})
    for k, v in b["diff"].items():
        print("  UQ -", k, [round(x, 4) for x in v["ci"]], "n_pos", v["n_pos"], {h: round(x, 3) for h, x in v["per_model"].items()})
    print("  cost m1", [round(x, 4) for x in b["cost_at_m1"]], "declared", {k: round(v, 4) for k, v in b["declared"].items()},
          "anchor09", round(b["anchor09"], 4))
    for n, c in b["curves"].items():
        print("  %-12s" % n, " ".join("%.3f" % x for x in c.values()))
    if "rules" in b:
        print("  rules", b["rules"], "core independent", b["core_generator_independent"])
