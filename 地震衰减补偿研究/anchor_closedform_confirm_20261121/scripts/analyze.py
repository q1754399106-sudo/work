"""Verdicts of PLAN_ZH.md section 3 (C1-C3, S1-S3) from results/*.json -> analysis/confirm.json and a printed summary."""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]; RES = HERE / "results"; OUT = HERE / "analysis"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import PROFILES, FAM_A  # noqa: E402

recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(RES.glob("*.json"))]
valid = [r for r in recs if r["valid"]]


def c1(rs):
    d = np.abs([r["F"]["pipe"]["alpha"] - r["kappa_hat"] for r in rs])
    return dict(n=len(d), median=float(np.median(d)), p90=float(np.percentile(d, 90)), frac_le_0_10=float(np.mean(d <= 0.10)),
                max=float(d.max()), pass_=bool(np.median(d) <= 0.03 and np.mean(d <= 0.10) >= 0.90))


prof = {}
for p in PROFILES:
    rs = [r for r in valid if r["profile"] == p]
    if not rs:
        prof[p] = dict(n=0); continue
    k = float(np.median([r["kappa_hat"] for r in rs])); g = lambda key: float(np.median([key(r) for r in rs]))
    F, U, M = g(lambda r: r["F"]["pipe"]["alpha"]), g(lambda r: r["U"]), g(lambda r: r["med"])
    prof[p] = dict(n=len(rs), family="A" if p in FAM_A else "B", kappa_hat=k, F=F, U=U, med=M, eF=F - k, eU=U - k, eM=M - k,
                   t_med=g(lambda r: r["t_med"]), share={b: g(lambda r, b=b: r["share_C"][b]) for b in ("tt", "tm", "mm")},
                   alpha_block={b: g(lambda r, b=b: r["alpha_block"][b]) for b in ("tt", "tm", "mm")},
                   khat_by_c={str(c): float(np.median([r["kappa_hat"] for r in rs if r["c_true"] == c])) for c in (0.7, 1.0, 1.4)
                              if any(r["c_true"] == c for r in rs)},
                   F_by_c={str(c): float(np.median([r["F"]["pipe"]["alpha"] for r in rs if r["c_true"] == c])) for c in (0.7, 1.0, 1.4)
                           if any(r["c_true"] == c for r in rs)})
pv = [v for v in prof.values() if v["n"]]
C1 = c1(valid); C3 = c1([r for r in valid if r["family"] == "B"]) if any(r["family"] == "B" for r in valid) else None
nC2 = sum(abs(v["eF"]) <= 0.05 for v in pv); nS2 = sum(abs(v["eM"]) <= 0.05 for v in pv)
nB = sum(r["family"] == "B" for r in recs); nBv = sum(r["family"] == "B" for r in valid)
undet = len(valid) < 0.8 * 102 or nBv < 0.8 * 36
out = dict(n_records=len(recs), n_valid=len(valid), n_B=nB, n_B_valid=nBv, undetermined=bool(undet),
           invalid=[f"{r['profile']}_{r['id']}" for r in recs if not r["valid"]],
           C1=C1, C2=dict(n_ok=int(nC2), n_profiles=len(pv), pass_=bool(nC2 >= 15)), C3=C3,
           PASS=bool(not undet and C1["pass_"] and nC2 >= 15 and C3 and C3["pass_"]),
           S1=dict(F_best=int(sum(abs(v["eF"]) < min(abs(v["eU"]), abs(v["eM"])) for v in pv)),
                   max_abs=dict(F=max(abs(v["eF"]) for v in pv), U=max(abs(v["eU"]) for v in pv), med=max(abs(v["eM"]) for v in pv))),
           S2=dict(n_ok=int(nS2), holds=bool(nS2 >= 15)),
           worst=sorted(([f"{r['profile']}_{r['id']}", r["F"]["pipe"]["alpha"] - r["kappa_hat"]] for r in valid),
                        key=lambda x: -abs(x[1]))[:10],
           profiles=prof)
OUT.mkdir(exist_ok=True); (OUT / "confirm.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"records {len(recs)} valid {len(valid)} (B {nBv}/{nB}) undetermined {undet} PASS {out['PASS']}")
print("C1", {k: round(v, 3) if isinstance(v, float) else v for k, v in C1.items()})
print("C3", None if C3 is None else {k: round(v, 3) if isinstance(v, float) else v for k, v in C3.items()})
print("C2", nC2, "/", len(pv), " S2", nS2, "/", len(pv), " S1 F best", out["S1"]["F_best"], out["S1"]["max_abs"])
print("profile fam n  khat    F     eF     eU     eM   t_med  share tt/tm/mm   alpha tt/tm/mm")
for p, v in prof.items():
    if v["n"]:
        print(f"{p:4s} {v['family']} {v['n']:2d} {v['kappa_hat']:.3f} {v['F']:.3f} {v['eF']:+.3f} {v['eU']:+.3f} {v['eM']:+.3f} {v['t_med']:.2f}"
              f"  {v['share']['tt']:+.2f}/{v['share']['tm']:+.2f}/{v['share']['mm']:+.2f}"
              f"  {v['alpha_block']['tt']:.2f}/{v['alpha_block']['tm']:.2f}/{v['alpha_block']['mm']:.2f}")
print("worst", [(a, round(b, 3)) for a, b in out["worst"]]); print("invalid", out["invalid"])
