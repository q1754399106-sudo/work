"""Rules R1-R4 and attributions A-C of PLAN_MMR_ZH.md -> analysis/mmr.json (and a printed summary)."""
import json, glob, collections
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results"
M_ALL = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0); UNDER = (1.4, 1.7, 2.0); NEAR = (0.85, 1.0, 1.2)
RULES = ("UQW", "MMJ", "RGJ", "RGC", "RGJ2"); rng_seed = 20261116


def load():
    semi, fd, info = collections.defaultdict(dict), collections.defaultdict(lambda: collections.defaultdict(list)), collections.defaultdict(list)
    for p in glob.glob(str(R / "mmr" / "*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8"))
        for r in J["records"]:
            info[(r["method"], "FD" if J["id"].startswith("FD") else "semi")].append(r["info"])
            if J["id"].startswith("FD"):
                fd[J["id"].split("_")[1]][(r["method"], r["m"])].append(r["S"])
            else:
                semi[J["id"]][(r["method"], r["m"])] = r["S"]
    for p in glob.glob(str(R / "stage0" / "*__VDB_m*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8"))
        if J["id"] in semi:
            semi[J["id"]][("VDB", J["m"])] = J["score"]["S"]
    for p in glob.glob(str(R / "stage2b" / "*_14dB_r*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8"))
        for r in J["records"]:
            if r["method"] == "VDB":
                fd[J["model"]][("VDB", r["m"])].append(r["S"])
    for p in glob.glob(str(R / "stageW" / "*.json")):
        J = json.loads(Path(p).read_text(encoding="utf-8")); fdj = J["id"].startswith("FD")
        if fdj and "14dB" not in J["id"]:
            continue
        for r in J["records"]:
            if r["method"] == "W0" or (r["method"] == "UQW_g0" and abs(r["sig"] - 0.447) < 1e-9):
                key = ("W0" if r["method"] == "W0" else "UQW_stageW", r["m"])
                if fdj:
                    fd[J["id"].split("_")[1]][key].append(r["S"])
                else:
                    semi[J["id"]][key] = r["S"]
    fdm = {h: {k: float(np.mean(v)) for k, v in d.items()} for h, d in fd.items()}
    return dict(semi), fdm, info


def mean_over(d, meth, ms):
    return float(np.mean([d[(meth, m)] for m in ms]))


def boot(x, rng):
    x = np.asarray(x); b = x[rng.integers(0, len(x), (2000, len(x)))].mean(1)
    return [float(x.mean()), float(np.quantile(b, 0.05)), float(np.quantile(b, 0.95))]


def diff(units, a, b, ms, rng=None):
    v = [mean_over(u, a, ms) - mean_over(u, b, ms) for u in units.values()]
    return boot(v, rng) if rng is not None else [float(np.mean(v)), int(sum(x > 0 for x in v)), len(v)]


if __name__ == "__main__":
    semi, fdm, info = load(); rng = np.random.default_rng(rng_seed); out = dict(n_semi=len(semi), n_fd=len(fdm))
    sb = lambda a, b, ms: diff(semi, a, b, ms, np.random.default_rng(rng_seed))  # noqa: E731
    out["check_uqw_vs_stageW_max_abs"] = max(abs(u[("UQW", m)] - u[("UQW_stageW", m)]) for u in semi.values() for m in M_ALL)
    out["R1"] = sb("RGJ", "UQW", UNDER); out["R2_semi"] = sb("RGJ", "UQW", NEAR); out["R2_fd"] = diff(fdm, "RGJ", "UQW", NEAR)
    out["R3"] = sb("RGJ", "UQW", M_ALL); out["R4_vdb"] = diff(fdm, "RGJ", "VDB", M_ALL); out["R4_uqw"] = diff(fdm, "RGJ", "UQW", M_ALL)
    out["pass"] = dict(R1=out["R1"][1] > 0, R2=out["R2_semi"][0] >= -0.01 and out["R2_fd"][0] >= -0.01,
                       R3=out["R3"][0] >= 0 and out["R3"][1] >= -0.01,
                       R4=out["R4_vdb"][0] > 0 and out["R4_vdb"][1] >= 5 and out["R4_uqw"][0] >= -0.01)
    for tag, (a, b) in dict(A=("RGJ", "MMJ"), B=("RGJ", "RGC"), C=("RGJ2", "RGJ")).items():
        out["attr_" + tag] = dict(semi=sb(a, b, M_ALL), semi_under=sb(a, b, UNDER), fd=diff(fdm, a, b, M_ALL))
    meths = RULES + ("VDB", "W0")
    out["per_m"] = {grp: {k: [float(np.mean([u[(k, m)] for u in d.values()])) for m in M_ALL] for k in meths}
                    for grp, d in (("semi", semi), ("fd14", fdm))}
    out["mean"] = {grp: {k: float(np.mean(v)) for k, v in pm.items()} for grp, pm in out["per_m"].items()}
    pis = {}
    for (meth, grp), L in info.items():
        if meth in ("MMJ", "RGJ"):
            P = np.array([i["pi"] for i in L]); pis[f"{meth}_{grp}"] = dict(pi_mean=P.mean(0).round(3).tolist(), gap_max=max(i["gap"] for i in L),
                                                                         w_u2=[float(P[:, -1].min()), float(P[:, -1].max())], w_u05=[float(P[:, 0].min()), float(P[:, 0].max())])
        elif meth == "RGC":
            P = np.array([i["pi_mean"] for i in L]); pis[f"RGC_{grp}"] = dict(pi_mean=P.mean(0).round(3).tolist(), gap_rel_p99_max=max(i["gap_rel_p99"] for i in L))
        elif meth == "RGJ2":
            P = np.array([i["pi"] for i in L]); u = np.array(L[0]["u"]); top = np.argsort(-P.mean(0))[:4]
            pis[f"RGJ2_{grp}"] = dict(top=[[u[i].tolist(), float(P.mean(0)[i])] for i in top], gap_max=max(i["gap"] for i in L))
    out["priors"] = pis
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "mmr.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k not in ("per_m", "priors")}, indent=1))
    for grp, pm in out["per_m"].items():
        print(grp); [print("  %-5s" % k, " ".join("%.3f" % x for x in v), " mean %.3f" % np.mean(v)) for k, v in pm.items()]
    print(json.dumps(pis, indent=1))
