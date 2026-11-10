"""Rules T1-T4 and attributions of PLAN_CB_ZH.md (amended) -> analysis/cb.json (and a printed summary).
UQW, RGJ, VDB and W0 come from analyze_mmr.load (results/mmr, stages 0, 2b and W); UQWt and CB from results/cb."""
import sys, json, glob, collections
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_mmr import load as load_mmr, M_ALL, UNDER, NEAR, boot, mean_over  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; SEED = 20261117; OVER = (0.5, 0.6, 0.7)
NEW = ("UQWt", "CB0.75", "CB0.5", "CB0.25")


def load():
    semi, fdm, _ = load_mmr(); fdraw = collections.defaultdict(lambda: collections.defaultdict(list)); infos = collections.defaultdict(list)
    for p in glob.glob(str(ROOT / "results" / "cb" / "*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8")); isfd = J["id"].startswith("FD")
        for r in J["records"]:
            if r["method"] == "info":
                infos["FD" if isfd else "semi"].append(r["info"]); continue
            if isfd:
                fdraw[J["id"].split("_")[1]][(r["method"], r["m"])].append(r["S"])
            else:
                semi[J["id"]][(r["method"], r["m"])] = r["S"]
    for h, d in fdraw.items():
        fdm[h].update({k: float(np.mean(v)) for k, v in d.items()})
    return semi, fdm, infos


def d_semi(semi, a, b, ms):
    return boot([mean_over(u, a, ms) - mean_over(u, b, ms) for u in semi.values()], np.random.default_rng(SEED))


def d_fd(fdm, a, b, ms):
    v = [mean_over(u, a, ms) - mean_over(u, b, ms) for u in fdm.values()]; return [float(np.mean(v)), int(sum(x > 0 for x in v)), len(v)]


if __name__ == "__main__":
    semi, fdm, infos = load(); out = dict(n_semi=len(semi), n_fd=len(fdm), n_semi_cb=sum(("CB0.5", 1.0) in u for u in semi.values()))
    out["T1"] = d_semi(semi, "CB0.5", "UQW", UNDER); out["T2_semi"] = d_semi(semi, "CB0.5", "UQW", NEAR); out["T2_fd"] = d_fd(fdm, "CB0.5", "UQW", NEAR)
    out["T3_uqw"] = d_fd(fdm, "CB0.5", "UQW", M_ALL); out["T3_vdb"] = d_fd(fdm, "CB0.5", "VDB", M_ALL); out["T4"] = d_semi(semi, "CB0.5", "UQW", M_ALL)
    out["pass"] = dict(T1=out["T1"][1] > 0, T2=out["T2_semi"][0] >= -0.01 and out["T2_fd"][0] >= -0.01,
                       T3=out["T3_uqw"][0] >= -0.01 and out["T3_vdb"][0] > 0 and out["T3_vdb"][1] >= 5, T4=out["T4"][0] >= 0 and out["T4"][1] >= -0.01)
    out["attr"] = {f"{a}-{b}": dict(semi=d_semi(semi, a, b, M_ALL), semi_under=d_semi(semi, a, b, UNDER), semi_near=d_semi(semi, a, b, NEAR),
                                     semi_over=d_semi(semi, a, b, OVER), fd=d_fd(fdm, a, b, M_ALL), fd_near=d_fd(fdm, a, b, NEAR))
                   for a, b in (("UQWt", "UQW"), ("CB0.5", "UQWt"), ("CB0.25", "UQW"), ("CB0.75", "UQW"), ("CB0.5", "RGJ"))}
    meths = ("UQW",) + NEW + ("RGJ", "VDB", "W0")
    out["per_m"] = {g: {k: [float(np.mean([u[(k, m)] for u in d.values()])) for m in M_ALL] for k in meths} for g, d in (("semi", semi), ("fd14", fdm))}
    out["mean"] = {g: {k: float(np.mean(v)) for k, v in pm.items()} for g, pm in out["per_m"].items()}
    sol = {}
    for g, L in infos.items():
        for t in ("0.75", "0.5", "0.25"):
            runs = [i["runs"][t] for i in L]; lam = np.array([r["lam"] for r in runs])
            sol[f"{g}_tau{t}"] = dict(viol_rel_max=max(r["viol_rel"] for r in runs), lam_mean=lam.mean(0).round(3).tolist(),
                                     frac_lam_only_ends=float(np.mean([(np.abs(l[1:-1]) < 1e-6).all() for l in lam])), lam_sum_median=float(np.median(lam.sum(1))))
        sol[f"{g}_R"] = dict(R_mm=float(np.mean([i["R_mm"] for i in L])), R_T=float(np.mean([i["R_T"] for i in L])))
    out["solver"] = sol
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "cb.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k not in ("per_m", "solver")}, indent=1))
    for g, pm in out["per_m"].items():
        print(g); [print("  %-7s" % k, " ".join("%.3f" % x for x in v), " mean %.3f" % np.mean(v)) for k, v in pm.items()]
    print(json.dumps(sol, indent=1))
