"""Analysis of PLAN_STAGEW_ZH.md: rules S1, S2, F1, F2, N1 per gamma variant; selection; descriptive. -> analysis/stageW.json"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import DEV, M_ALL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; RW = ROOT / "results" / "stageW"
MODELS = [f"H{i}" for i in range(1, 7)]; RR = (1, 2, 3); SIG = 0.447
rng = np.random.default_rng(20261114); IDX = rng.integers(0, len(DEV), (2000, len(DEV)))
W = {d: json.loads((RW / f"{d}.json").read_text(encoding="utf-8"))["records"] for d in DEV}
WF = {(h, lv, r): json.loads((RW / f"FD_{h}_{lv}_r{r}.json").read_text(encoding="utf-8"))["records"]
      for h in MODELS for lv in ("14dB", "08dB") for r in RR}
S1R = {d: json.loads((ROOT / "results" / "stage1" / f"{d}.json").read_text(encoding="utf-8"))["records"] for d in DEV}
S2B = {k: json.loads((ROOT / "results" / "stage2b" / f"{k[0]}_{k[1]}_r{k[2]}.json").read_text(encoding="utf-8"))["records"] for k in WF}


def pick(recs, method, m, sig=None, gamma=None):
    hit = [r["S"] for r in recs if r["method"] == method and r["m"] == m and (sig is None or abs(r["sig"] - sig) < 1e-12)
           and (gamma is None or r.get("gamma") == gamma)]
    assert len(hit) == 1, (method, m, sig, gamma, hit); return hit[0]


def s0(d, method, m):
    return json.loads((ROOT / "results" / "stage0" / f"{d}__{method}_m{m:g}.json").read_text())["score"]["S"]


def semi_mean(get):                      # per dataset, mean over the 9 m
    return np.array([np.mean([get(d, m) for m in M_ALL]) for d in DEV])


def fd_mean(get, lv):                    # per model, mean over m and realisations
    return np.array([np.mean([get((h, lv, r), m) for m in M_ALL for r in RR]) for h in MODELS])


def ci(x):
    b = x[IDX].mean(axis=1); lo, hi = np.quantile(b, [0.05, 0.95]); return [float(x.mean()), float(lo), float(hi)]


semi = dict(UQIQ=semi_mean(lambda d, m: pick(S1R[d], "UQ_lognormal", m, SIG)), VDB=semi_mean(lambda d, m: s0(d, "VDB", m)),
            IQ70=semi_mean(lambda d, m: s0(d, "IQ70", m)), C1=semi_mean(lambda d, m: pick(S1R[d], "C1", m)),
            W0=semi_mean(lambda d, m: pick(W[d], "W0", m)))
fdm = {lv: dict(UQIQ=fd_mean(lambda k, m: pick(S2B[k], "UQ", m, SIG), lv), VDB=fd_mean(lambda k, m: pick(S2B[k], "VDB", m), lv),
                IQ70=fd_mean(lambda k, m: pick(S2B[k], "IQ70", m), lv), C1=fd_mean(lambda k, m: pick(S2B[k], "C1", m), lv),
                W0=fd_mean(lambda k, m: pick(WF[k], "W0", m), lv)) for lv in ("14dB", "08dB")}
out = dict(variants={})
for g in (0.0, 1.0):
    U = semi_mean(lambda d, m: pick(W[d], "UQW_g%d" % g, m, SIG)); F = {lv: fd_mean(lambda k, m: pick(WF[k], "UQW_g%d" % g, m, SIG), lv)
                                                                        for lv in ("14dB", "08dB")}
    v = dict(semi_mean=float(U.mean()), fd14_mean=float(F["14dB"].mean()), fd08_mean=float(F["08dB"].mean()),
             S1=float((U - semi["UQIQ"]).mean()), S2=ci(U - semi["VDB"]),
             F1=float((F["14dB"] - fdm["14dB"]["VDB"]).mean()),
             F2=dict(mean=float((F["14dB"] - fdm["14dB"]["C1"]).mean()), n_pos=int(((F["14dB"] - fdm["14dB"]["C1"]) > 0).sum())),
             N1_semi=ci(U - semi["W0"]),
             N1_fd=dict(mean=float((F["14dB"] - fdm["14dB"]["W0"]).mean()), n_pos=int(((F["14dB"] - fdm["14dB"]["W0"]) > 0).sum())),
             desc=dict(semi_minus_IQ70=float((U - semi["IQ70"]).mean()), fd14_minus_UQIQ=float((F["14dB"] - fdm["14dB"]["UQIQ"]).mean()),
                       fd08_minus_VDB=float((F["08dB"] - fdm["08dB"]["VDB"]).mean()),
                       fd14_per_model_minus_VDB=dict(zip(MODELS, map(float, F["14dB"] - fdm["14dB"]["VDB"])))))
    v["pass"] = dict(S1=v["S1"] >= -0.01, S2=v["S2"][1] > 0, F1=v["F1"] >= -0.01, F2=v["F2"]["mean"] > 0 and v["F2"]["n_pos"] >= 5,
                     N1=v["N1_semi"][1] > 0 and v["N1_fd"]["mean"] > 0 and v["N1_fd"]["n_pos"] >= 5)
    v["all"] = all(v["pass"].values()); out["variants"]["gamma%d" % g] = v
ok = [k for k, v in out["variants"].items() if v["all"]]
out["selected"] = max(ok, key=lambda k: out["variants"][k]["semi_mean"] + out["variants"][k]["fd14_mean"]) if ok else None
out["baselines"] = dict(semi={k: float(v.mean()) for k, v in semi.items()}, fd={lv: {k: float(v.mean()) for k, v in d.items()} for lv, d in fdm.items()})
cur = {}
for name, meth, sig in (("W0", "W0", 0.0), ("UQW_g0", "UQW_g0", SIG), ("UQW_g1", "UQW_g1", SIG), ("UQW_g0_0.245", "UQW_g0", 0.245),
                        ("UQW_g0_0.894", "UQW_g0", 0.894)):
    cur[name] = dict(semi={f"{m:g}": float(np.mean([pick(W[d], meth, m, sig) for d in DEV])) for m in M_ALL},
                     fd14={f"{m:g}": float(np.mean([pick(WF[(h, "14dB", r)], meth, m, sig) for h in MODELS for r in RR])) for m in M_ALL})
out["curves"] = cur
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "stageW.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("baselines", json.dumps(out["baselines"], indent=None))
for k, v in out["variants"].items():
    print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a not in ("desc",)})
    print("   desc", v["desc"])
print("selected", out["selected"])
for n, c in cur.items():
    print("%-14s semi" % n, " ".join("%.3f" % x for x in c["semi"].values()), "| fd14", " ".join("%.3f" % x for x in c["fd14"].values()))
