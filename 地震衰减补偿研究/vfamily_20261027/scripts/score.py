"""Scoring of PLAN_ZH.md ("指标", "判定", "只描述"). One row = scene x noise x Q error.
SC / QScan + SC from cache/solved; Standard, QScan + Standard, FullCalRefit, QScanRefit from
cache/runner/results. Relative waveform error ||X-T||/||T||. Judged per class (W layered, X complex).
Writes analysis/rows.csv and analysis/decisions.json.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
sys.path.insert(0, str(RESEARCH / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(ROOT / "cache" / "inputs" / "novelty_evidence_20261025" / "scripts"))
import common as cmn  # noqa: E402
import csv, hashlib, json  # noqa: E402,E401
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

HASHES = json.loads((ROOT / "logs" / "MODEL_HASHES.json").read_text(encoding="utf-8"))
CLASSES = {"W": [f"W{k}" for k in range(1, 7)], "X": [f"X{k}" for k in range(1, 4)]}
NOISES = ["Moderate", "Strong"]
ERRORS = [-30, 0, 30]
ARCHIVED = ["Standard", "QScan_Standard", "FullCalRefit", "QScanRefit"]
METHODS = ["QScan_SC", "SC", *ARCHIVED]

rows, infos = [], []
for cls, scenes in CLASSES.items():
    for m in scenes:
        tp = ROOT / "truth" / f"{m}_truth.mat"
        if f"truth/{m}_truth.mat" not in HASHES:
            print(f"{m}: not generated, skipped")
            continue
        assert hashlib.sha256(tp.read_bytes()).hexdigest() == HASHES[f"truth/{m}_truth.mat"], m
        tr = loadmat(tp, simplify_cells=True)["tr"]
        for noise in NOISES:
            T = np.asarray(next(t for t in tr if t["name"] == noise)["target"], float)
            nT = max(np.linalg.norm(T), 1e-30)
            for e in ERRORS:
                res = loadmat(ROOT / "cache/runner/results" / f"{m}_Q{e:+d}_{noise}.mat", simplify_cells=True)
                out = {k: np.asarray(res["gathers"][k], float) for k in ARCHIVED}
                row = dict(cls=cls, scene=m, noise=noise, q_error=e, qscan_factor=float(res["scan"]["factor"]))
                for label, key in (("SC", "SC"), ("QSC", "QScan_SC")):
                    f = ROOT / "cache/solved" / f"{m}_{noise}_Q{e:+d}_{label}"
                    out[key] = np.load(f.with_suffix(".npz"))["out"]
                    info = json.loads(f.with_suffix(".json").read_text(encoding="utf-8"))
                    infos.append(info)
                    row[f"converged_{key}"] = bool(info["converged"])
                for k in METHODS:
                    row[f"rel_{k}"] = float(np.linalg.norm(out[k] - T) / nT)
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
    side = "SC side" if ab["stable"] else "Standard side" if ba["stable"] else "none"
    return dict(stable_side=side, sc_side=ab, standard_side=ba)


sel = lambda cls, wrong: [x for x in rows if x["cls"] == cls and (x["q_error"] != 0) == wrong]
conv = lambda rr, key: [x for x in rr if x[f"converged_{key}"]]
decisions, conclusions = {}, {}
for cls in CLASSES:
    if not sel(cls, True):
        continue
    P1 = both(sel(cls, True), "QScan_SC", "QScan_Standard")
    P2 = both(sel(cls, False), "SC", "Standard")
    sides = (P1["stable_side"], P2["stable_side"])
    conclusions[cls] = ("reproduced the reversal" if sides == ("Standard side", "Standard side")
                        else "did not reproduce the reversal" if "SC side" in sides else "partial")
    decisions[cls] = dict(P1_Q_wrong=P1, P2_Q_correct=P2, conclusion=conclusions[cls],
                          sensitivity_converged_only=dict(P1=both(conv(sel(cls, True), "QScan_SC"), "QScan_SC", "QScan_Standard"),
                                                          P2=both(conv(sel(cls, False), "SC"), "SC", "Standard")))
vals = list(conclusions.values())
overall = ("supports the scene-family explanation" if "reproduced the reversal" in vals
           else "refutes the scene-family explanation" if vals and all(v == "did not reproduce the reversal" for v in vals)
           else "inconclusive")
per_scene = {}
for m in sorted({x["scene"] for x in rows}):
    w_, r_ = [x for x in rows if x["scene"] == m and x["q_error"]], [x for x in rows if x["scene"] == m and not x["q_error"]]
    per_scene[m] = dict(QScanSC_over_QScanStandard=float(np.exp(np.mean([np.log(x["rel_QScan_SC"] / x["rel_QScan_Standard"]) for x in w_]))),
                        SC_over_Standard=float(np.exp(np.mean([np.log(x["rel_SC"] / x["rel_Standard"]) for x in r_]))),
                        QScanSC_over_QScanRefit=float(np.exp(np.mean([np.log(x["rel_QScan_SC"] / x["rel_QScanRefit"]) for x in w_]))),
                        qscan_factors=sorted({x["qscan_factor"] for x in rows if x["scene"] == m}))
cmn.writejson(an / "decisions.json", dict(
    classes=decisions, overall=overall, solves=len(infos),
    unconverged=[dict(scene=i["model"], noise=i["noise"], q_error=i["q_error"], label=i["label"], kkt=i["kkt_ratio"])
                 for i in infos if not i["converged"]],
    descriptive=dict(per_scene=per_scene,
                     class_means={c: {f"{k}_{lab}": float(np.mean([x[f"rel_{k}"] for x in sel(c, lab == "wrong")]))
                                      for k in METHODS for lab in ("wrong", "right")} for c in decisions})))
for m, v in per_scene.items():
    print(f"{m}: QScan+SC/QScan+Std {v['QScanSC_over_QScanStandard']:.3f}  SC/Std {v['SC_over_Standard']:.3f}")
for c, d in decisions.items():
    p1, p2 = d["P1_Q_wrong"], d["P2_Q_correct"]
    print(f"{c}: P1 {p1['stable_side']} (ratio {p1['sc_side']['mean_ratio']:.3f}, SC side wins {p1['sc_side']['wins']}/{p1['sc_side']['n']}) | "
          f"P2 {p2['stable_side']} (ratio {p2['sc_side']['mean_ratio']:.3f}, {p2['sc_side']['wins']}/{p2['sc_side']['n']}) -> {d['conclusion']}")
print("overall:", overall)
