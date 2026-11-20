"""Analysis of the development run (PLAN_ZH.md sections 3-5) -> analysis/dev.json and a printed summary."""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; QC = ROOT.parent / "qunknown_comp_20261109" / "results"
M_ALL = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0)
SEG = dict(over=(0.5, 0.6, 0.7), near=(0.85, 1.0, 1.2), under=(1.4, 1.7, 2.0), all=M_ALL)
RULES = ("SIR", "SIU", "UQW", "UQWl", "CBL0.5", "RGJ", "VDB", "W0")
rng = np.random.default_rng(20261120)


def load():
    semi, fdr, info, chk = {}, {}, {}, {}
    for p in sorted((ROOT / "results" / "dev").glob("*.json")):
        r = json.loads(p.read_text(encoding="utf-8")); recs = [x for x in r["records"] if x["method"] != "info"]
        d = {(x["method"], x["m"]): x for x in recs}
        info[r["id"]] = {x["m"]: x["info"] for x in r["records"] if x["method"] == "info"}
        if r["id"].startswith("FD_"):
            fdr[r["id"]] = d
        else:
            semi[r["id"]] = d
        old = {}
        for sub in ("mmr", "cbl"):
            q = QC / sub / p.name
            if q.exists():
                for x in json.loads(q.read_text(encoding="utf-8"))["records"]:
                    if x["method"] in ("UQW", "RGJ", "CBL0.5", "UQWl"):
                        old[(x["method"], x["m"])] = x["S"]
        for k, v in old.items():
            if k in d:
                chk[k[0]] = max(chk.get(k[0], 0.0), abs(d[k]["S"] - v))
    return semi, fdr, info, chk


def table(units, metric):
    """{rule: (n_units, 9) array}"""
    return {r: np.array([[units[u][(r, m)][metric] for m in M_ALL] for u in units]) for r in RULES}


def seg(a, s):
    idx = [M_ALL.index(m) for m in SEG[s]]; return a[:, idx].mean(1)


def boot(v):
    b = rng.integers(0, len(v), (2000, len(v))); bm = v[b].mean(1); return [float(v.mean()), float(np.quantile(bm, 0.05)), float(np.quantile(bm, 0.95))]


def fd_models(fdr):
    """average the 3 realisations per model -> pseudo units keyed by model"""
    out = {}
    for h in sorted({k.split("_")[1] for k in fdr}):
        ks = [k for k in fdr if k.split("_")[1] == h]
        out[h] = {key: {mt: float(np.mean([fdr[k][key][mt] for k in ks])) for mt in ("S", "M1", "M1_lag")} for key in fdr[ks[0]]}
    return out


if __name__ == "__main__":
    semi, fdr, info, chk = load(); fdm = fd_models(fdr)
    res = dict(n_semi=len(semi), n_fd=len(fdr), check_max_abs_diff=chk)
    TS, TT, TF = table(semi, "S"), table(semi, "T2"), table(fdm, "S")
    res["curves"] = {k: {r: np.round(t[r].mean(0), 4).tolist() for r in RULES} for k, t in
                     (("semi_S", TS), ("semi_T2", TT), ("fd_S", TF), ("semi_M1", table(semi, "M1")), ("semi_M1lag", table(semi, "M1_lag")),
                      ("fd_M1", table(fdm, "M1")), ("fd_M1lag", table(fdm, "M1_lag")))}
    res["trueQ_W0_m1"] = dict(semi_S=float(TS["W0"][:, 4].mean()), semi_T2=float(TT["W0"][:, 4].mean()), fd_S=float(TF["W0"][:, 4].mean()))
    diffs = {}
    for X in ("SIR", "SIU"):
        for Y in ("CBL0.5", "UQW", "RGJ", "VDB", "W0"):
            for s in SEG:
                diffs[f"{X}-{Y}|{s}"] = dict(semi_S=boot(seg(TS[X] - TS[Y], s)), semi_T2=boot(seg(TT[X] - TT[Y], s)),
                                             fd_S=[float(seg(TF[X] - TF[Y], s).mean()), int((seg(TF[X] - TF[Y], s) > 0).sum())])
    for s in SEG:
        diffs[f"CBL0.5-UQW|{s}"] = dict(semi_S=boot(seg(TS["CBL0.5"] - TS["UQW"], s)), semi_T2=boot(seg(TT["CBL0.5"] - TT["UQW"], s)),
                                         fd_S=[float(seg(TF["CBL0.5"] - TF["UQW"], s).mean()), int((seg(TF["CBL0.5"] - TF["UQW"], s) > 0).sum())])
    res["diffs"] = diffs
    gates = {}
    for X in ("SIR", "SIU"):
        d = lambda Y, s: diffs[f"{X}-{Y}|{s}"]
        g = dict(D1=d("CBL0.5", "all")["semi_S"][0] >= 0.010 and d("CBL0.5", "all")["semi_S"][1] > 0,
                 D2=d("CBL0.5", "all")["fd_S"][0] >= -0.010,
                 D3=d("UQW", "near")["semi_S"][0] >= -0.010 and d("UQW", "near")["fd_S"][0] >= -0.010,
                 D4=d("CBL0.5", "all")["semi_T2"][0] >= -0.010,
                 D5=d("VDB", "all")["fd_S"][0] > 0 and d("VDB", "all")["fd_S"][1] >= 5)
        g["pass"] = all(g.values()); gates[X] = g
    res["gates"] = gates
    sir = [(u, m, v["SIR"]) for u, mm in info.items() for m, v in mm.items()]
    res["SIR_solver"] = dict(gap_max=max(x[2]["gap"] for x in sir), its_max=max(x[2]["its"] for x in sir),
                             sec_mean=float(np.mean([x[2]["sec"] for x in sir])),
                             pi_semi_mean=np.round(np.mean([x[2]["pi"] for x in sir if not x[0].startswith("FD")], 0), 3).tolist(),
                             pi_fd_mean=np.round(np.mean([x[2]["pi"] for x in sir if x[0].startswith("FD")], 0), 3).tolist())
    res["alpha_median"] = {r: dict(semi=np.round(np.median([np.array(v[r]["alpha_q"])[1] for u, mm in info.items() if not u.startswith("FD") for v in mm.values()], 0), 3).tolist(),
                                   fd=np.round(np.median([np.array(v[r]["alpha_q"])[1] for u, mm in info.items() if u.startswith("FD") for v in mm.values()], 0), 3).tolist())
                           for r in ("SIU", "SIR")}
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "dev.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print("units semi %d fd %d; check vs stored S, max |diff| per rule:" % (len(semi), len(fdr)), chk)
    for k in ("semi_S", "semi_T2", "fd_S"):
        print(k); [print("  %-7s" % r, " ".join("%.3f" % v for v in res["curves"][k][r]), " mean %.3f" % np.mean(res["curves"][k][r])) for r in RULES]
    for k, v in diffs.items():
        if k.split("|")[0] in ("SIR-CBL0.5", "SIU-CBL0.5", "SIR-UQW", "SIU-UQW", "CBL0.5-UQW", "SIR-RGJ", "SIR-VDB", "SIU-VDB"):
            print("%-18s semi S %+.4f [%+.4f,%+.4f]  T2 %+.4f [%+.4f,%+.4f]  FD S %+.4f (%d/6)" % (k, *v["semi_S"], *v["semi_T2"], *v["fd_S"]))
    print("gates", gates); print("SIR solver", res["SIR_solver"]); print("alpha", res["alpha_median"]); print("trueQ", res["trueQ_W0_m1"])
