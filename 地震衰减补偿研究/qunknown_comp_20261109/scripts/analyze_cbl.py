"""Rules T1-T4 and attributions of PLAN_CBL_ZH.md -> analysis/cbl.json (and a printed summary).
UQW, RGJ, VDB, W0 (results/mmr, stages 0, 2b, W) and UQWt, CB0.5 (results/cb) via analyze_cb.load; UQWl and CBL from results/cbl.
results/cb must be present in the working tree (restore it from git if the sparse checkout hid it)."""
import sys, json, glob, collections
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_cb import load as load_cb, d_semi, d_fd, OVER  # noqa: E402
from analyze_mmr import M_ALL, UNDER, NEAR  # noqa: E402
from cbl import UL, CON  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; NEW = ("UQWl", "CBL0.5", "CBL0.25")


def load():
    semi, fdm, _ = load_cb(); fdraw = collections.defaultdict(lambda: collections.defaultdict(list)); infos = collections.defaultdict(list)
    for p in glob.glob(str(ROOT / "results" / "cbl" / "*.json")):
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


if __name__ == "__main__":
    semi, fdm, infos = load()
    out = dict(n_semi=len(semi), n_fd=len(fdm), n_semi_cbl=sum(("CBL0.5", 1.0) in u for u in semi.values()),
               n_semi_cb=sum(("CB0.5", 1.0) in u for u in semi.values()))
    out["T1"] = d_semi(semi, "CBL0.5", "UQW", UNDER); out["T2_semi"] = d_semi(semi, "CBL0.5", "UQW", NEAR); out["T2_fd"] = d_fd(fdm, "CBL0.5", "UQW", NEAR)
    out["T3_uqw"] = d_fd(fdm, "CBL0.5", "UQW", M_ALL); out["T3_vdb"] = d_fd(fdm, "CBL0.5", "VDB", M_ALL); out["T4"] = d_semi(semi, "CBL0.5", "UQW", M_ALL)
    out["pass"] = dict(T1=out["T1"][1] > 0, T2=out["T2_semi"][0] >= -0.01 and out["T2_fd"][0] >= -0.01,
                       T3=out["T3_uqw"][0] >= -0.01 and out["T3_vdb"][0] > 0 and out["T3_vdb"][1] >= 5, T4=out["T4"][0] >= 0 and out["T4"][1] >= -0.01)
    out["attr"] = {f"{a}-{b}": dict(semi=d_semi(semi, a, b, M_ALL), semi_under=d_semi(semi, a, b, UNDER), semi_near=d_semi(semi, a, b, NEAR),
                                     semi_over=d_semi(semi, a, b, OVER), fd=d_fd(fdm, a, b, M_ALL), fd_near=d_fd(fdm, a, b, NEAR))
                   for a, b in (("UQWl", "UQW"), ("CBL0.5", "UQWl"), ("CBL0.5", "CB0.5"), ("CBL0.25", "UQW"), ("UQWl", "UQWt"), ("CBL0.5", "RGJ"))}
    meths = ("UQW", "UQWt", "CB0.5") + NEW + ("RGJ", "VDB", "W0")
    out["per_m"] = {g: {k: [float(np.mean([u[(k, m)] for u in d.values()])) for m in M_ALL] for k in meths} for g, d in (("semi", semi), ("fd14", fdm))}
    out["mean"] = {g: {k: float(np.mean(v)) for k, v in pm.items()} for g, pm in out["per_m"].items()}
    sol = {}; ucon = UL[CON]
    for g, L in infos.items():
        for t in ("0.5", "0.25"):
            runs = [i["runs"][t] for i in L]; lam = np.array([r["lam"] for r in runs])
            sol[f"{g}_tau{t}"] = dict(viol_rel_max=max(r["viol_rel"] for r in runs), lam_mean=lam.mean(0).round(3).tolist(), u=ucon.round(3).tolist(),
                                     frac_lam_only_ends=float(np.mean([(np.abs(l[1:-1]) < 1e-6).all() for l in lam])), lam_sum_median=float(np.median(lam.sum(1))),
                                     its_max=max(r["its"] for r in runs))
        sol[f"{g}_R"] = dict(R_mm=float(np.mean([i["R_mm"] for i in L])), R_T=float(np.mean([i["R_T"] for i in L])))
    out["solver"] = sol
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "cbl.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k not in ("per_m", "solver")}, indent=1))
    for g, pm in out["per_m"].items():
        print(g); [print("  %-7s" % k, " ".join("%.3f" % x for x in v), " mean %.3f" % np.mean(v)) for k, v in pm.items()]
    print(json.dumps(sol, indent=1))
