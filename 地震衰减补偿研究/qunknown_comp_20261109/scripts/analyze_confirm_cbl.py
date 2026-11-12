"""Analysis of PLAN_CONFIRM_CBL_ZH.md: rules C1-C6 and descriptive tables -> analysis/confirm_cbl.json.
UQW, VDB, IQ70, W0 from results/confirmW; UQWl, CBL0.5, CBL0.25, UQW_re from results/confirm_cbl.
Semi-synthetic unit: dataset (12), bootstrap 2000, 90 %, seed 20261118; FD unit: model (6), mean over 3 realisations."""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage0 import M_ALL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]; RW = ROOT / "results" / "confirmW"; RC = ROOT / "results" / "confirm_cbl"
IDS5 = [f"Z{k}_c{c}_r{r}" for k in "MNO" for c in (0.75, 1.05) for r in (0, 1)]
MODELS2 = [f"H{i}" for i in range(7, 13)]; R = (1, 2, 3)
OLD = ("UQW", "VDB", "IQ70", "W0"); NEW = ("UQWl", "CBL0.5", "CBL0.25", "UQW_re")
SEG = dict(all=M_ALL, under=(1.4, 1.7, 2.0), near=(0.85, 1.0, 1.2), over=(0.5, 0.6, 0.7))
STATES = (("semi", "cal"), ("semi", "nd"), ("fd", "30dB"), ("fd", "14dB"))
IDX = np.random.default_rng(20261118).integers(0, len(IDS5), (2000, len(IDS5)))
_cache = {}


def recs(path):
    if path not in _cache:
        _cache[path] = {(x["method"], x["m"]): x["S"] for x in json.loads(path.read_text(encoding="utf-8"))["records"] if x["method"] != "info"}
    return _cache[path]


def S(path_old, path_new, method, m):
    return recs(path_old if method in OLD else path_new)[(method, m)]


def units(state, method, ms):
    kind, v = state
    if kind == "semi":
        return np.array([np.mean([S(RW / f"S5_{d}_{v}.json", RC / f"S5_{d}_{v}.json", method, m) for m in ms]) for d in IDS5])
    return np.array([np.mean([S(RW / f"FD2_{h}_{v}_r{r}.json", RC / f"FD2_{h}_{v}_r{r}.json", method, m) for m in ms for r in R]) for h in MODELS2])


def summ(state, a, b, ms):
    x = units(state, a, ms) - units(state, b, ms); r = dict(mean=float(x.mean()), npos=int((x > 0).sum()), n=len(x))
    if state[0] == "semi":
        r["lo"], r["hi"] = [float(q) for q in np.quantile(x[IDX].mean(1), [0.05, 0.95])]
    return r


if __name__ == "__main__":
    st = lambda s: f"{s[0]}_{s[1]}"; out = dict(diffs={}, per_m={}, checks={})
    for s in STATES:
        out["diffs"][st(s)] = {f"{a}-{b}_{g}": summ(s, a, b, ms) for a, b in (("CBL0.5", "UQW"), ("CBL0.5", "VDB"), ("UQWl", "UQW"), ("CBL0.5", "UQWl"),
                                                                           ("CBL0.25", "UQW")) for g, ms in SEG.items()}
        out["per_m"][st(s)] = {k: [float(units(s, k, (m,)).mean()) for m in M_ALL] for k in OLD + NEW}
        out["checks"][st(s)] = float(max(np.abs(units(s, "UQW_re", (m,)) - units(s, "UQW", (m,))).max() for m in M_ALL))
    D = out["diffs"]; d = lambda s, key: D[st(s)][key]
    semi_cal, semi_nd, fd30, fd14 = STATES; C = {}
    C["C1"] = d(semi_cal, "CBL0.5-UQW_under")["lo"] > 0
    x = d(fd14, "CBL0.5-UQW_under"); C["C2"] = x["mean"] > 0 and x["npos"] >= 5
    C["C3"] = all((d(s, "CBL0.5-UQW_all")["lo"] if s[0] == "semi" else d(s, "CBL0.5-UQW_all")["mean"]) >= -0.01 for s in STATES)
    C["C4"] = all(d(s, "CBL0.5-UQW_near")["mean"] >= -0.01 for s in STATES)
    C["C5"] = all(d(s, "CBL0.5-UQW_over")["mean"] >= -0.02 for s in STATES)
    x = d(fd14, "CBL0.5-VDB_all"); C["C6"] = d(semi_nd, "CBL0.5-VDB_all")["lo"] > 0 and x["mean"] > 0 and x["npos"] >= 5
    out["pass"] = C; out["confirmed"] = all(C.values())
    out["mean"] = {k: {m: float(np.mean(v)) for m, v in pm.items()} for k, pm in out["per_m"].items()}
    (ROOT / "analysis").mkdir(exist_ok=True); (ROOT / "analysis" / "confirm_cbl.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(dict(pass_=C, confirmed=out["confirmed"], checks=out["checks"]), indent=1))
    for s in STATES:
        print(st(s))
        for k, v in D[st(s)].items():
            print("  %-22s" % k, " ".join("%s=%.4f" % (a, b) if isinstance(b, float) else "%s=%s" % (a, b) for a, b in v.items()))
        for k, v in out["per_m"][st(s)].items():
            print("  %-8s" % k, " ".join("%.3f" % x for x in v), " mean %.3f" % np.mean(v))
