"""Verdicts of experiment 12 (PLAN_ZH.md section 4) from results/synth/*.json; writes analysis/summary.json."""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
from scipy import stats  # noqa: E402
recs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((ROOT / "results" / "synth").glob("T_*.json"))]
assert len(recs) == 20, len(recs)
PROP = "SRSCa_1.3+ph"


def paired(a, b, key):
    d = np.array([r["scores"][a][key] - r["scores"][b][key] for r in recs]); m = float(d.mean())
    h = float(stats.t.ppf(0.95, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d)))
    return dict(mean=m, lo=m - h, hi=m + h, positive=int((d > 0).sum()), n=len(d))


cmp = {}
for a, b in ((PROP, "B2_SCa"), ("SCaK", "B2_SCa"), (PROP, "B1_IQ"), (PROP, "B3_VDB"), ("IQd", "B1_IQ"),
             (PROP, "SRSCa"), (PROP, "SRSCa_1.3"), (PROP, "SCa_true"), ("SCaK", "B1_IQ"), ("SCaK", "B3_VDB"),
             ("B1_IQ", "B3_VDB")):
    if a in recs[0]["scores"] and b in recs[0]["scores"]:
        cmp[f"{a} - {b}"] = {k: paired(a, b, k) for k in ("S", "M1", "M1_lag")}
v = {
    "H1": cmp[f"{PROP} - B2_SCa"]["S"]["mean"] > 0 and cmp[f"{PROP} - B2_SCa"]["S"]["lo"] > 0,
    "H1b": "SCaK - B2_SCa" in cmp and cmp["SCaK - B2_SCa"]["S"]["mean"] > 0 and cmp["SCaK - B2_SCa"]["S"]["lo"] > 0,
    "H2": cmp[f"{PROP} - B1_IQ"]["S"]["lo"] > -0.02,
    "H3": cmp[f"{PROP} - B3_VDB"]["S"]["lo"] > -0.02,
    "H4": cmp["IQd - B1_IQ"]["S"]["mean"] > 0 and cmp["IQd - B1_IQ"]["S"]["lo"] > 0,
}
names = list(recs[0]["scores"])
mean = {m: {k: float(np.mean([r["scores"][m][k] for r in recs])) for k in ("M1", "M1_lag", "S")} for m in names}
byc = {}
for c in sorted({r["c_true"] for r in recs}):
    rr = [r for r in recs if r["c_true"] == c]
    byc[str(c)] = {m: float(np.mean([r["scores"][m]["S"] for r in rr])) for m in names}
csr = np.array([r["sr"]["c"] / r["c_true"] for r in recs])
sr = dict(abs_log_err=float(np.mean(np.abs(np.log(csr)))), mean_ratio=float(csr.mean()),
          within25=float(np.mean(np.abs(csr - 1) <= 0.25)))
out = dict(verdicts=v, comparisons=cmp, mean=mean, by_c=byc, sr=sr, lines=recs[0]["lines"])
(ROOT / "analysis" / "summary.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps(v)); print(json.dumps(sr))
for m in sorted(names, key=lambda m: -mean[m]["S"]):
    print(f"{m:16s} M1 {mean[m]['M1']:.3f}  M1_lag {mean[m]['M1_lag']:.3f}  S {mean[m]['S']:.3f}")
for k, x in cmp.items():
    print(k, {kk: f"{y['mean']:+.3f} [{y['lo']:+.3f}, {y['hi']:+.3f}] {y['positive']}/{y['n']}" for kk, y in x.items()})
