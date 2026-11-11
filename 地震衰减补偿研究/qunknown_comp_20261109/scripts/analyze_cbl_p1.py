"""Summary of the P1 development check (cbl_p1.py) -> analysis/cbl_p1.json: CBL0.5 - UQW, - VDB, UQWl - UQW, CBL0.5 - UQWl
by segment on semi DEV 4x noise (24 datasets, bootstrap 2000, 90 %, seed 20261118) and FD H1-H6 30 dB (model means over 3 realisations)."""
import json, sys, glob, collections
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import M_ALL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SEG = dict(all=M_ALL, under=(1.4, 1.7, 2.0), near=(0.85, 1.0, 1.2), over=(0.5, 0.6, 0.7))
IDX = np.random.default_rng(20261118).integers(0, 24, (2000, 24))

if __name__ == "__main__":
    semi = collections.defaultdict(dict); fd = collections.defaultdict(lambda: collections.defaultdict(list)); sol = collections.defaultdict(list)
    for p in sorted(glob.glob(str(ROOT / "results" / "cbl_p1" / "*.json"))):
        J = json.loads(Path(p).read_text(encoding="utf-8")); isfd = J["id"].startswith("FD")
        for r in J["records"]:
            if r["method"] == "info":
                sol["fd" if isfd else "semi"].append(max(v["viol_rel"] for v in r["info"]["runs"].values()))
            elif isfd:
                fd[J["id"].split("_")[1]][(r["method"], r["m"])].append(r["S"])
            else:
                semi[J["id"]][(r["method"], r["m"])] = r["S"]
    fdm = {h: {k: float(np.mean(v)) for k, v in d.items()} for h, d in fd.items()}; out = dict(n_semi=len(semi), n_fd=len(fdm))
    mo = lambda u, k, ms: np.mean([u[(k, m)] for m in ms])
    for g, D in (("semi_x4", semi), ("fd_30dB", fdm)):
        res = {}
        for a, b in (("CBL0.5", "UQW"), ("CBL0.5", "VDB"), ("UQWl", "UQW"), ("CBL0.5", "UQWl")):
            for s, ms in SEG.items():
                x = np.array([mo(u, a, ms) - mo(u, b, ms) for u in D.values()]); r = dict(mean=float(x.mean()), npos=int((x > 0).sum()), n=len(x))
                if g == "semi_x4" and len(x) == 24:
                    r["lo"], r["hi"] = [float(q) for q in np.quantile(x[IDX].mean(1), [0.05, 0.95])]
                res[f"{a}-{b}_{s}"] = r
        res["per_m"] = {k: [float(np.mean([u[(k, m)] for u in D.values()])) for m in M_ALL] for k in ("UQW", "UQWl", "CBL0.5", "CBL0.25", "VDB")}
        out[g] = res
    out["viol_rel_max"] = {g: float(max(v)) for g, v in sol.items()}
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "cbl_p1.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for g in ("semi_x4", "fd_30dB"):
        print(g)
        for k, v in out[g].items():
            if k != "per_m":
                print("  %-22s" % k, " ".join("%s=%.4f" % (a, b) if isinstance(b, float) else "%s=%s" % (a, b) for a, b in v.items()))
        for k, v in out[g]["per_m"].items():
            print("  %-8s" % k, " ".join("%.3f" % x for x in v), " mean %.3f" % np.mean(v))
    print("viol_rel_max", out["viol_rel_max"])
