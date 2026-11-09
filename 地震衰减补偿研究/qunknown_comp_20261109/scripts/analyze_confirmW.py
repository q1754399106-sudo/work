"""Analysis of PLAN_CONFIRM_W_ZH.md: rules K1-K6; descriptive tables. -> analysis/confirmW.json"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import M_ALL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; RC = ROOT / "results" / "confirmW"
IDS5 = [f"Z{k}_c{c}_r{r}" for k in "MNO" for c in (0.75, 1.05) for r in (0, 1)]
MODELS2 = [f"H{i}" for i in range(7, 13)]; R = (1, 2, 3)
METHODS = ("UQW", "W0", "VDB", "IQ70", "UQIQ", "C1")
rng = np.random.default_rng(20261116); IDX = rng.integers(0, len(IDS5), (2000, len(IDS5)))


def get(path, method, m):
    hit = [x["S"] for x in json.loads(path.read_text(encoding="utf-8"))["records"] if x["method"] == method and x["m"] == m]
    assert len(hit) == 1, (path, method, m); return hit[0]


def semi_units(version, method):
    return np.array([np.mean([get(RC / f"S5_{d}_{version}.json", method, m) for m in M_ALL]) for d in IDS5])


def fd_units(lv, method):
    return np.array([np.mean([get(RC / f"FD2_{h}_{lv}_r{r}.json", method, m) for m in M_ALL for r in R]) for h in MODELS2])


def ci(x):
    b = x[IDX].mean(axis=1); lo, hi = np.quantile(b, [0.05, 0.95]); return [float(x.mean()), float(lo), float(hi)]


U = {("semi", v): {k: semi_units(v, k) for k in METHODS} for v in ("cal", "nd")}
U.update({("fd", lv): {k: fd_units(lv, k) for k in METHODS} for lv in ("30dB", "14dB")})
diff = lambda s, a, b: U[s][a] - U[s][b]
semi_pos = lambda x: ci(x)[1] > 0
fd_pos = lambda x: x.mean() > 0 and (x > 0).sum() >= 5
out = dict(mean_S={f"{a}_{b}": {k: float(v.mean()) for k, v in d.items()} for (a, b), d in U.items()}, diffs={})
for (a, b) in U:
    out["diffs"][f"{a}_{b}"] = {f"UQW-{k}": (ci(diff((a, b), "UQW", k)) if a == "semi" else
                                             dict(mean=float(diff((a, b), "UQW", k).mean()), n_pos=int((diff((a, b), "UQW", k) > 0).sum()),
                                                  per_model=dict(zip(MODELS2, map(float, diff((a, b), "UQW", k))))))
                                for k in METHODS if k != "UQW"}
K = dict(K1=float(diff(("semi", "cal"), "UQW", "VDB").mean()) >= -0.01,
         K2=float(diff(("fd", "30dB"), "UQW", "VDB").mean()) >= -0.01,
         K3=semi_pos(diff(("semi", "nd"), "UQW", "VDB")),
         K4=fd_pos(diff(("fd", "14dB"), "UQW", "VDB")))
for kk, ref in (("K5", "IQ70"), ("K6", "W0")):
    K[kk] = all([semi_pos(diff(("semi", v), "UQW", ref)) for v in ("cal", "nd")] + [fd_pos(diff(("fd", lv), "UQW", ref)) for lv in ("30dB", "14dB")])
out["rules"] = {k: bool(v) for k, v in K.items()}; out["confirmed"] = all(K.values())
out["curves"] = {f"{a}_{b}": {k: {f"{m:g}": float(np.mean([get(RC / (f"S5_{d}_{b}.json" if a == "semi" else f"FD2_{d}_{b}_r{r}.json"), k, m)
                                                            for d in (IDS5 if a == "semi" else MODELS2) for r in ((0,) if a == "semi" else R)]))
                                  for m in M_ALL} for k in METHODS} for (a, b) in U}
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "confirmW.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for s, v in out["mean_S"].items():
    print(s, {k: round(x, 4) for k, x in v.items()})
for s, v in out["diffs"].items():
    print(s, {k: ([round(x, 4) for x in d] if isinstance(d, list) else dict(mean=round(d["mean"], 4), n_pos=d["n_pos"])) for k, d in v.items()})
print("rules", out["rules"], "confirmed", out["confirmed"])
for s, c in out["curves"].items():
    for k in ("UQW", "VDB", "W0", "IQ70"):
        print("%-10s %-5s" % (s, k), " ".join("%.3f" % x for x in c[k].values()))
