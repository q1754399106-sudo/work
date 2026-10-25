"""Scoring of PLAN_ZH.md ("指标", "判定"). Relative waveform error ||X-T||/||T|| (primary) and
normalized squared error per replicate; scene values are means over r = 1..6.
QScan+SC from cache/solved; SC (nominal Q) from sc_neutral_20261026/cache/solved; Standard and
QScan+Standard from the fd2d sweep cache. Every local input is checked against cache/SNAPSHOT.json.
Writes analysis/replicates.csv, analysis/scenes.csv and analysis/decisions.json.
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

LOCAL = json.loads((ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))["local_sha256"]
MODELS = [f"H{m}" for m in range(1, 7)]
CARBONATE = {"H2", "H4", "H6"}
SCENES = ["Q100_14dB", "Q100_08dB", "Q070_14dB"]
METHODS = ["QScan_SC", "SC", "QScan_Standard", "Standard"]


def checked(rel):
    p = RESEARCH / rel
    assert hashlib.sha256(p.read_bytes()).hexdigest() == LOCAL[rel], f"input changed: {rel}"
    return p


reps, infos = [], []
for m in MODELS:
    T = np.asarray(loadmat(ROOT / "cache/inputs/fd2d_20261025/truth" / f"{m}_truth.mat")["target"], float)
    nT = np.linalg.norm(T)
    for s in SCENES:
        for r in range(7):
            L = loadmat(checked(f"fd2d_20261025/cache/sweep/{m}_{s}/r{r:03d}.mat"), simplify_cells=True)
            best = int(L["info"]["best"])
            out = {"Standard": L["STD"][:, :, 0], "QScan_Standard": L["STD"][:, :, best],
                   "SC": np.load(checked(f"sc_neutral_20261026/cache/solved/{m}_{s}_r{r}_SC.npz"))["out"]}
            f = ROOT / "cache/solved" / f"{m}_{s}_r{r}_QSC"
            out["QScan_SC"] = np.load(f.with_suffix(".npz"))["out"]
            info = json.loads(f.with_suffix(".json").read_text())
            assert info["qscan_best"] == best
            infos.append(dict(model=m, scene=s, r=r, factor=info["qscan_factor"], converged=info["converged"],
                              kkt=info["kkt_ratio"], iterations=info["iterations"]))
            for meth in METHODS:
                e = np.linalg.norm(np.asarray(out[meth], float) - T) / nT
                reps.append(dict(model=m, scene=s, r=r, method=meth, rel=e, nmse=e * e))

an = ROOT / "analysis"
an.mkdir(exist_ok=True)
with open(an / "replicates.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(reps[0]))
    w.writeheader()
    w.writerows(reps)


def val(m, s, meth, key="rel", rs=range(1, 7)):
    return float(np.mean([x[key] for x in reps if x["model"] == m and x["scene"] == s and x["method"] == meth and x["r"] in rs]))


scenes = []
for m in MODELS:
    for s in SCENES:
        row = dict(model=m, scene=s, carbonate=m in CARBONATE)
        for meth in METHODS:
            row[f"rel_{meth}"] = val(m, s, meth)
            row[f"nmse_{meth}"] = val(m, s, meth, "nmse")
            row[f"rel0_{meth}"] = val(m, s, meth, rs=[0])
        row["qsc_converged_r1_6"] = all(i["converged"] for i in infos if i["model"] == m and i["scene"] == s and i["r"] >= 1)
        scenes.append(row)
with open(an / "scenes.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(scenes[0]))
    w.writeheader()
    w.writerows(scenes)


def gate(rows, a, b):
    """Is method a stably better than method b (PLAN_ZH.md 判定)?"""
    ma, mb = np.mean([x[f"rel_{a}"] for x in rows]), np.mean([x[f"rel_{b}"] for x in rows])
    ratios = [x[f"rel_{a}"] / x[f"rel_{b}"] for x in rows]
    wins = sum(x[f"rel_{a}"] < x[f"rel_{b}"] for x in rows)
    c1, c2, c3 = ma <= 0.95 * mb, wins >= int(np.ceil(0.75 * len(rows))), max(ratios) <= 1.05
    return dict(mean_a=ma, mean_b=mb, mean_ratio=ma / mb, wins=wins, n=len(rows), worst_ratio=max(ratios),
                c1_mean=bool(c1), c2_wins=bool(c2), c3_worst=bool(c3), stable=bool(c1 and c2 and c3))


def both(rows, a, b):
    ab, ba = gate(rows, a, b), gate(rows, b, a)
    verdict = f"{a} stably better than {b}" if ab["stable"] else f"{b} stably better than {a}" if ba["stable"] else "no stable advantage"
    return dict(verdict=verdict, a_better=ab, b_better=ba)


wrong = [x for x in scenes if x["scene"] == "Q070_14dB"]
right = [x for x in scenes if x["scene"] != "Q070_14dB"]
H1 = both(wrong, "QScan_SC", "QScan_Standard")
H2 = both(right, "SC", "QScan_SC")
conv = lambda rows: [x for x in rows if x["qsc_converged_r1_6"]]
desc = dict(
    Q070_QScanSC_vs_SC=both(wrong, "QScan_SC", "SC"),
    Q0_QScanSC_vs_QScanStandard=both(right, "QScan_SC", "QScan_Standard"),
    Q0_14dB_SC_vs_QScanSC=both([x for x in right if x["scene"] == "Q100_14dB"], "SC", "QScan_SC"),
    Q0_08dB_SC_vs_QScanSC=both([x for x in right if x["scene"] == "Q100_08dB"], "SC", "QScan_SC"),
    Q0_carbonate_SC_vs_QScanSC=both([x for x in right if x["carbonate"]], "SC", "QScan_SC"),
    Q0_no_carbonate_SC_vs_QScanSC=both([x for x in right if not x["carbonate"]], "SC", "QScan_SC"),
    mean_rel={s: {k: float(np.mean([x[f"rel_{k}"] for x in scenes if x["scene"] == s])) for k in METHODS} for s in SCENES},
    mean_rel_r0={s: {k: float(np.mean([x[f"rel0_{k}"] for x in scenes if x["scene"] == s])) for k in METHODS} for s in SCENES},
)
cmn.writejson(an / "decisions.json", dict(
    H1_Q_wrong=H1, H2_Q_correct=H2, qsc_solves=len(infos),
    qsc_unconverged=[i for i in infos if not i["converged"]],
    sensitivity_converged_only=dict(H1=both(conv(wrong), "QScan_SC", "QScan_Standard") if conv(wrong) else None,
                                    n_H1=len(conv(wrong)),
                                    H2=both(conv(right), "SC", "QScan_SC") if conv(right) else None,
                                    n_H2=len(conv(right))),
    descriptive=desc))
for x in scenes:
    print(f"{x['model']} {x['scene']}: QScan+SC {x['rel_QScan_SC']:.4f} SC {x['rel_SC']:.4f} "
          f"QScan+Std {x['rel_QScan_Standard']:.4f} Std {x['rel_Standard']:.4f} conv {x['qsc_converged_r1_6']}")
print("H1:", H1["verdict"], "| H2:", H2["verdict"])
