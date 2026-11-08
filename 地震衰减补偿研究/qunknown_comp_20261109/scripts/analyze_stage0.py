"""Analysis of PLAN_STAGE0_ZH.md: mean loss L(m) = S(1) - S(m) per method, rules R1-R3 on the better method M*.
-> analysis/stage0.json"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import DEV, M_ALL, M_SC  # noqa: E402

R = ROOT / "results" / "stage0"; MS = {"IQ70": M_ALL, "VDB": M_ALL, "SCa": M_SC}
S = {meth: {m: np.array([json.loads((R / f"{d}__{meth}_m{m:g}.json").read_text())["score"]["S"] for d in DEV]) for m in ms}
     for meth, ms in MS.items()}
rng = np.random.default_rng(20261109); NB = 2000; IDX = rng.integers(0, len(DEV), (NB, len(DEV)))


def ci(x):
    b = x[IDX].mean(axis=1); lo, hi = np.quantile(b, [0.05, 0.95]); return float(x.mean()), float(lo), float(hi)


out = dict(mean_S={}, loss={})
for meth, ms in MS.items():
    out["mean_S"][meth] = {f"{m:g}": float(S[meth][m].mean()) for m in ms}
    out["loss"][meth] = {f"{m:g}": ci(S[meth][1.0] - S[meth][m]) for m in ms if m != 1.0}
best = max(("IQ70", "SCa"), key=lambda k: S[k][1.0].mean()); L = lambda m: S[best][1.0] - S[best][m]
sym = ci(0.5 * (L(0.7) + L(1.4))); asym = ci(L(0.7) - L(1.4)); far = ci(0.5 * (L(0.5) + L(2.0))); fasym = ci(L(0.5) - L(2.0))
R1 = sym[0] >= 0.03 or asym[0] >= 0.03; R2 = (not R1) and far[0] >= 0.03
out["rules"] = dict(M_star=best, sym_07_14=sym, asym_07_minus_14=asym, far_05_20=far, asym_05_minus_20=fasym,
                    R1=R1, R2=R2, R3=not (R1 or R2))
# descriptive: VDB against the better tuned method at each m, and the actual SR error of experiment 17
out["vdb_minus_best"] = {f"{m:g}": ci(S["VDB"][m] - S[best][m]) for m in M_SC}
dev = {d: json.loads((BASE / "mu_tune_20261031" / "results" / "dev" / f"{d}.json").read_text()) for d in DEV}
ratio = np.array([dev[d]["c_sr"] / dev[d]["c_eff_cum"] for d in DEV])
out["sr_ratio"] = dict(median=float(np.median(ratio)), min=float(ratio.min()), max=float(ratio.max()),
                       abs_log_median=float(np.median(np.abs(np.log(ratio)))))
# descriptive: S(m = 1, this run) minus the SR-route score of experiment 17 (environment difference ~1.5e-4)
out["oracle_minus_sr"] = {k: ci(S[k][1.0] - np.array([dev[d]["scores"][ref]["S"] for d in DEV]))
                          for k, ref in (("IQ70", "IQ70"), ("SCa", "SR_mu0.25"))}
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "stage0.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for meth, ms in MS.items():
    print(meth, " ".join(f"m{m:g}:{out['mean_S'][meth][f'{m:g}']:.3f}" for m in ms))
for meth in MS:
    print("loss", meth, " ".join(f"m{k}:{v[0]:+.3f}[{v[1]:+.3f},{v[2]:+.3f}]" for k, v in out["loss"][meth].items()))
print("VDB - best", {k: round(v[0], 3) for k, v in out["vdb_minus_best"].items()})
print("SR ratio", out["sr_ratio"]); print("oracle - SR", out["oracle_minus_sr"]); print("rules", out["rules"])
