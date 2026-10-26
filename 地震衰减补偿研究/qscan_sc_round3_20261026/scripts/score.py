"""Scoring of PLAN_ZH.md ("指标", "判定", "只描述"). One row = gather x Q error (one noise realization).
SC / QScan + SC from cache/solved; Standard, QScan + Standard, FullCalRefit and QScanRefit from the
archived round3 results. Relative waveform error ||X-T||/||T|| (primary) and normalized squared error.
Writes analysis/rows.csv and analysis/decisions.json.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
SNAP = ROOT / "cache" / "inputs"
sys.path.insert(0, str(RESEARCH / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(SNAP / "novelty_evidence_20261025" / "scripts"))
import common as cmn  # noqa: E402
import csv, hashlib, json, re  # noqa: E402,E401
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

S = json.loads((ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))
GROUPS = {"Kennett": ["F3", "F4"], "V": ["V4", "V5"]}
NOISES = ["Moderate", "Strong"]
ERRORS = [-50, -30, -10, 0, 10, 30, 50]
ARCHIVED = {"Standard": "Standard", "QScan_Standard": "QScan_Standard", "FullCalRefit": "FullCalRefit",
            "QScanRefit": "QScanRefit"}
METHODS = ["QScan_SC", "SC", *ARCHIVED]


def snap(rel):
    p = SNAP / rel
    assert hashlib.sha256(p.read_bytes()).hexdigest() == S["files"][rel]["sha256"], f"input changed: {rel}"
    return loadmat(p, simplify_cells=True)


rows, infos = [], []
for group, sets in GROUPS.items():
    for name in sets:
        tr = snap(f"round3_20261024/dev/{name}_truth.mat")["tr"]
        tr = tr if isinstance(tr, list) else [tr]
        for noise in NOISES:
            T = np.asarray(next(t for t in tr if t["name"] == noise)["target"], float)
            nT = max(np.linalg.norm(T), 1e-30)
            for e in ERRORS:
                g = snap(f"round3_20261024/results/{name}_Q{e:+d}_{noise}.mat")["gathers"]
                out = {k: np.asarray(g[v], float) for k, v in ARCHIVED.items()}
                row = dict(group=group, set=name, noise=noise, q_error=e)
                for label, key in (("SC", "SC"), ("QSC", "QScan_SC")):
                    f = ROOT / "cache/solved" / f"{name}_{noise}_Q{e:+d}_{label}"
                    out[key] = np.load(f.with_suffix(".npz"))["out"]
                    info = json.loads(f.with_suffix(".json").read_text(encoding="utf-8"))
                    infos.append(info)
                    row[f"converged_{key}"] = bool(info["converged"])
                    row[f"seconds_{key}"] = float(info["seconds"])
                    row["qscan_factor"] = float(info["qscan_factor"])
                for k in METHODS:
                    err = float(np.linalg.norm(out[k] - T) / nT)
                    row[f"rel_{k}"] = err
                    row[f"nmse_{k}"] = err * err
                rows.append(row)

an = ROOT / "analysis"
an.mkdir(exist_ok=True)
with open(an / "rows.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)


def gate(rr, a, b):
    """Is method a stably better than method b (PLAN_ZH.md 判定)?"""
    ma, mb = np.mean([x[f"rel_{a}"] for x in rr]), np.mean([x[f"rel_{b}"] for x in rr])
    ratios = [x[f"rel_{a}"] / x[f"rel_{b}"] for x in rr]
    wins = sum(x[f"rel_{a}"] < x[f"rel_{b}"] for x in rr)
    c1, c2, c3 = ma <= 0.95 * mb, wins >= int(np.ceil(0.75 * len(rr))), max(ratios) <= 1.05
    return dict(mean_a=float(ma), mean_b=float(mb), mean_ratio=float(ma / mb), wins=int(wins), n=len(rr),
                worst_ratio=float(max(ratios)), c1_mean=bool(c1), c2_wins=bool(c2), c3_worst=bool(c3),
                stable=bool(c1 and c2 and c3))


def both(rr, a, b):
    ab, ba = gate(rr, a, b), gate(rr, b, a)
    verdict = f"{a} stably better than {b}" if ab["stable"] else f"{b} stably better than {a}" if ba["stable"] else "no stable advantage"
    return dict(verdict=verdict, a_better=ab, b_better=ba)


sel = lambda **kw: [x for x in rows if all((x[k] in v) if isinstance(v, (list, tuple, set)) else x[k] == v for k, v in kw.items())]
wrong = [e for e in ERRORS if e]
conv = lambda rr: [x for x in rr if x["converged_SC"] and x["converged_QScan_SC"]]
R1 = {g: both(sel(group=g, q_error=wrong), "QScan_SC", "QScan_Standard") for g in GROUPS}
k_ok, v_ok = R1["Kennett"]["a_better"]["stable"], R1["V"]["a_better"]["stable"]
k_rev, v_rev = R1["Kennett"]["b_better"]["stable"], R1["V"]["b_better"]["stable"]
R1_conclusion = ("conflict with H1-H6: QScan_Standard stably better on " + ",".join(g for g, r in (("Kennett", k_rev), ("V", v_rev)) if r)
                 if k_rev or v_rev else "confirmed on 1-D generators" if k_ok and v_ok
                 else "generator-dependent: confirmed on Kennett only" if k_ok else "not confirmed")
R2 = both(sel(q_error=0), "SC", "QScan_SC")

desc = dict(
    by_group_error={g: {str(e): {k: float(np.mean([x[f"rel_{k}"] for x in sel(group=g, q_error=e)])) for k in METHODS}
                        for e in ERRORS} for g in GROUPS},
    by_set_wrong={s: both(sel(set=s, q_error=wrong), "QScan_SC", "QScan_Standard") for g in GROUPS.values() for s in g},
    rescue_QScanSC_vs_SC={g: both(sel(group=g, q_error=wrong), "QScan_SC", "SC") for g in GROUPS},
    QScanSC_vs_QScanRefit={g: both(sel(group=g, q_error=wrong), "QScan_SC", "QScanRefit") for g in GROUPS},
    Q0_QScanSC_vs_QScanStandard=both(sel(q_error=0), "QScan_SC", "QScan_Standard"),
    Q0_by_group_SC_vs_QScanSC={g: both(sel(group=g, q_error=0), "SC", "QScan_SC") for g in GROUPS},
    all_rows_QScanSC_vs_QScanStandard=both(rows, "QScan_SC", "QScan_Standard"),
)
# External plan, stage-1 pass criterion 1: "QScan+Standard mean error better than or equal to QScan+SC".
desc["external_criterion_error"] = {g: dict(mean_QScan_Standard=float(np.mean([x["rel_QScan_Standard"] for x in sel(group=g)])),
                                            mean_QScan_SC=float(np.mean([x["rel_QScan_SC"] for x in sel(group=g)])))
                                    for g in GROUPS}
for g, d in desc["external_criterion_error"].items():
    d["met"] = d["mean_QScan_Standard"] <= d["mean_QScan_SC"]

# Runtime (descriptive only): SC seconds per solve vs MATLAB seconds per round3 row (QScan scan,
# two refits and coherence shrinkage, 4 threads), from the cumulative times in the round3 logs.
matlab = {}
for log in ("combo_F34.log", "combo_V4F.log", "combo_V5.log"):
    prev = 0.0
    for line in (SNAP / "round3_20261024/logs" / log).read_text(encoding="utf-8-sig").splitlines():
        m = re.match(r"(\w+) Q([+-]\d+) (\w+) factor [\d.]+ done ([\d.]+)s", line)
        if m:
            matlab[(m[1], int(m[2]), m[3])] = float(m[4]) - prev
            prev = float(m[4])
desc["runtime_seconds"] = {s: dict(SC_per_solve=float(np.mean([x["seconds_SC"] for x in sel(set=s)])),
                                   QScanSC_per_solve=float(np.mean([x["seconds_QScan_SC"] for x in sel(set=s)])),
                                   matlab_row_all_routes=float(np.mean([matlab[(s, e, n)] for e in ERRORS for n in NOISES])))
                           for g in GROUPS.values() for s in g}

cmn.writejson(an / "decisions.json", dict(
    R1_Q_wrong=R1, R1_conclusion=R1_conclusion, R2_Q_correct=R2, solves=len(infos),
    unconverged=[dict(set=i["set"], noise=i["noise"], q_error=i["q_error"], label=i["label"], kkt=i["kkt_ratio"])
                 for i in infos if not i["converged"]],
    novelty_consistency={f"{i['set']}_{i['noise']}": i["novelty_relative_difference"] for i in infos
                         if "novelty_relative_difference" in i},
    sensitivity_converged_only=dict(
        R1={g: both(conv(sel(group=g, q_error=wrong)), "QScan_SC", "QScan_Standard") for g in GROUPS},
        R2=both(conv(sel(q_error=0)), "SC", "QScan_SC") if conv(sel(q_error=0)) else None),
    descriptive=desc))
for x in rows:
    print(f"{x['set']} {x['noise']:8s} Q{x['q_error']:+3d}: QScan+SC {x['rel_QScan_SC']:.4f} SC {x['rel_SC']:.4f} "
          f"QScan+Std {x['rel_QScan_Standard']:.4f} Std {x['rel_Standard']:.4f}")
for g, r in R1.items():
    print(f"R1-{g}:", r["verdict"], f"ratio {r['a_better']['mean_ratio']:.3f} wins {r['a_better']['wins']}/{r['a_better']['n']}")
print("R1:", R1_conclusion, "| R2:", R2["verdict"])
