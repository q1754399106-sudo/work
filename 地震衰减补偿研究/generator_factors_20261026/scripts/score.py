"""Scoring of PLAN_ZH.md ("指标", "判定", "只描述").
Per row (model x noise x Q error): relative waveform error of SC / QScan + SC (cache/solved) and of
Standard, QScan + Standard, FullCalRefit, QScanRefit (cache/runner/results). Per model:
L_wrong = mean over the 4 Q-wrong rows of ln(QScan+SC / QScan+Standard); L_0 = mean over the 2 Q-right
rows of ln(SC / Standard). Writes analysis/rows.csv, analysis/models.csv and analysis/decisions.json.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
sys.path.insert(0, str(RESEARCH / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(ROOT / "cache" / "inputs" / "novelty_evidence_20261025" / "scripts"))
import common as cmn  # noqa: E402
import csv, hashlib, itertools, json  # noqa: E402,E401
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

HASHES = json.loads((ROOT / "logs" / "MODEL_HASHES.json").read_text(encoding="utf-8"))
SEEDS = [1, 2, 3]
CELLS = list(itertools.product([0, 1], repeat=3))
NOISES = ["Moderate", "Strong"]
ERRORS = [-30, 0, 30]
ARCHIVED = ["Standard", "QScan_Standard", "FullCalRefit", "QScanRefit"]
METHODS = ["QScan_SC", "SC", *ARCHIVED]
FACTORS = ["A_sparse", "B_structured", "C_solver_operator"]


def name(s, c):
    return f"Z{s}_{c[0]}{c[1]}{c[2]}"


rows, models, infos = [], [], []
for s in SEEDS:
    for c in CELLS:
        m = name(s, c)
        tp = ROOT / "truth" / f"{m}_truth.mat"
        assert hashlib.sha256(tp.read_bytes()).hexdigest() == HASHES[f"truth/{m}_truth.mat"], m
        tr = loadmat(tp, simplify_cells=True)["tr"]
        mr = []
        for noise in NOISES:
            T = np.asarray(next(t for t in tr if t["name"] == noise)["target"], float)
            nT = max(np.linalg.norm(T), 1e-30)
            for e in ERRORS:
                res = loadmat(ROOT / "cache/runner/results" / f"{m}_Q{e:+d}_{noise}.mat", simplify_cells=True)
                out = {k: np.asarray(res["gathers"][k], float) for k in ARCHIVED}
                row = dict(model=m, seed=s, A_sparse=c[0], B_structured=c[1], C_solver_operator=c[2], noise=noise,
                           q_error=e, qscan_factor=float(res["scan"]["factor"]))
                for label, key in (("SC", "SC"), ("QSC", "QScan_SC")):
                    f = ROOT / "cache/solved" / f"{m}_{noise}_Q{e:+d}_{label}"
                    out[key] = np.load(f.with_suffix(".npz"))["out"]
                    info = json.loads(f.with_suffix(".json").read_text(encoding="utf-8"))
                    infos.append(info)
                    row[f"converged_{key}"] = bool(info["converged"])
                for k in METHODS:
                    row[f"rel_{k}"] = float(np.linalg.norm(out[k] - T) / nT)
                rows.append(row)
                mr.append(row)
        wrong = [x for x in mr if x["q_error"]]
        right = [x for x in mr if not x["q_error"]]
        models.append(dict(model=m, seed=s, A_sparse=c[0], B_structured=c[1], C_solver_operator=c[2],
                           L_wrong=float(np.mean([np.log(x["rel_QScan_SC"] / x["rel_QScan_Standard"]) for x in wrong])),
                           L_0=float(np.mean([np.log(x["rel_SC"] / x["rel_Standard"]) for x in right])),
                           converged_wrong=all(x["converged_QScan_SC"] for x in wrong),
                           converged_0=all(x["converged_SC"] for x in right)))

an = ROOT / "analysis"
an.mkdir(exist_ok=True)
for fn, data in (("rows.csv", rows), ("models.csv", models)):
    with open(an / fn, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0]))
        w.writeheader()
        w.writerows(data)

L = {(x["seed"], x["A_sparse"], x["B_structured"], x["C_solver_operator"]): x for x in models}


def cell_mean(c, key, keep=None):
    v = [L[(s, *c)][key] for s in SEEDS if keep is None or keep(L[(s, *c)])]
    return float(np.mean(v)) if v else None


def factor(i, key, keep=None):
    d = []
    for s in SEEDS:
        for other in itertools.product([0, 1], repeat=2):
            hi, lo = list(other), list(other)
            hi.insert(i, 1)
            lo.insert(i, 0)
            a, b = L[(s, *hi)], L[(s, *lo)]
            if keep is None or (keep(a) and keep(b)):
                d.append(a[key] - b[key])
    d = np.asarray(d)
    mean, pos, neg, n = float(d.mean()), int((d > 0).sum()), int((d < 0).sum()), len(d)
    need10, need9 = int(np.ceil(10 / 12 * n)), int(np.ceil(9 / 12 * n))
    if mean >= np.log(1.10) and pos >= need10:
        verdict = "main driver"
    elif mean >= np.log(1.05) and pos >= need9:
        verdict = "has effect"
    elif mean <= -np.log(1.05) and neg >= need9:
        verdict = "reverse effect"
    else:
        verdict = "no clear effect"
    return dict(verdict=verdict, mean_d=mean, ratio=float(np.exp(mean)), positive=pos, negative=neg, n=n,
                pairs=[float(x) for x in d])


def interaction(i, j, key):
    """Mean pair difference of factor i at each level of factor j (descriptive)."""
    outd = {}
    for lev in (0, 1):
        d = []
        for s in SEEDS:
            for c in CELLS:
                if c[i] == 1 and c[j] == lev:
                    lo = list(c)
                    lo[i] = 0
                    d.append(L[(s, *c)][key] - L[(s, *lo)][key])
        outd[f"{FACTORS[j]}={lev}"] = float(np.mean(d))
    return outd


rep_v, rep_f = cell_mean((1, 1, 1), "L_wrong"), cell_mean((0, 0, 0), "L_wrong")
reproduced = rep_v > 0 and rep_f < 0
cw = lambda x: x["converged_wrong"]
c0 = lambda x: x["converged_0"]
cmn.writejson(an / "decisions.json", dict(
    reproduction=dict(V4_like_111=rep_v, F_like_000=rep_f, reproduced=bool(reproduced)),
    primary_L_wrong={f: factor(i, "L_wrong") for i, f in enumerate(FACTORS)},
    secondary_L_0={f: factor(i, "L_0") for i, f in enumerate(FACTORS)},
    sensitivity_converged_only=dict(primary={f: factor(i, "L_wrong", cw) for i, f in enumerate(FACTORS)},
                                    secondary={f: factor(i, "L_0", c0) for i, f in enumerate(FACTORS)}),
    unconverged=[dict(model=i["model"], noise=i["noise"], q_error=i["q_error"], label=i["label"], kkt=i["kkt_ratio"])
                 for i in infos if not i["converged"]],
    descriptive=dict(
        cells={"".join(map(str, c)): dict(L_wrong=cell_mean(c, "L_wrong"), L_0=cell_mean(c, "L_0"),
                                           **{f"rel_{k}_{lab}": float(np.mean([x[f"rel_{k}"] for x in rows
                                                if (x["A_sparse"], x["B_structured"], x["C_solver_operator"]) == c
                                                and ((x["q_error"] != 0) == (lab == "wrong"))]))
                                              for k in METHODS for lab in ("wrong", "right")})
               for c in CELLS},
        interactions_L_wrong={f"{FACTORS[i]}|{FACTORS[j]}": interaction(i, j, "L_wrong")
                              for i in range(3) for j in range(3) if i != j},
        interactions_L_0={f"{FACTORS[i]}|{FACTORS[j]}": interaction(i, j, "L_0")
                          for i in range(3) for j in range(3) if i != j},
        qscan_factors={m["model"]: sorted({x["qscan_factor"] for x in rows if x["model"] == m["model"]}) for m in models})))
for x in models:
    print(f"{x['model']}: L_wrong {x['L_wrong']:+.3f} (x{np.exp(x['L_wrong']):.3f})  L_0 {x['L_0']:+.3f} (x{np.exp(x['L_0']):.3f})")
print(f"reproduction: 111 {rep_v:+.3f}, 000 {rep_f:+.3f} -> {reproduced}")
for i, f in enumerate(FACTORS):
    p, q = factor(i, "L_wrong"), factor(i, "L_0")
    print(f"{f}: L_wrong {p['verdict']} (x{p['ratio']:.3f}, {p['positive']}/{p['n']} >0) | "
          f"L_0 {q['verdict']} (x{q['ratio']:.3f}, {q['positive']}/{q['n']} >0)")
