"""Scoring of PLAN_ZH.md ("指标", "判定"). Relative waveform error ||X-T||/||T|| (primary) and
normalized squared error ||X-T||^2/||T||^2 over all 99 traces and 451 samples, per replicate; scene
values are means over r = 1..6. Standard, FullCalRefit and QScan+Standard come from the fd2d sweep
cache (factor index 0 = nominal Q; QScan's choice = 1-based grid index 'best'), SC and SC0 from
cache/solved. Writes analysis/replicates.csv, analysis/scenes.csv and analysis/decisions.json.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
sys.path.insert(0, str(ROOT / "cache" / "inputs" / "novelty_evidence_20261025" / "scripts"))
import common as cmn  # noqa: E402
import csv, json  # noqa: E402,E401
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

# The repository reorganisation c625fd6 moved fd2d under 地震衰减补偿研究/ (PLAN_ADDENDUM_ZH.md);
# every sweep file is checked against the SHA-256 recorded in cache/SNAPSHOT.json before use.
SWEEP = ROOT.parent / "地震衰减补偿研究" / "fd2d_20261025" / "cache" / "sweep"
_snap = json.loads((ROOT / "cache" / "SNAPSHOT.json").read_text(encoding="utf-8"))["fd2d_sweep_sha256"]
import hashlib  # noqa: E402
for _k, _h in _snap.items():
    _p = SWEEP / Path(_k).relative_to("fd2d_20261025/cache/sweep")
    assert hashlib.sha256(_p.read_bytes()).hexdigest() == _h, f"sweep file changed: {_p}"
print(f"{len(_snap)} fd2d sweep files verified against the pre-run snapshot", flush=True)
MODELS = [f"H{m}" for m in range(1, 7)]
CARBONATE = {"H2", "H4", "H6"}
SCENES = ["Q100_14dB", "Q100_08dB", "Q070_14dB"]
PRIMARY = ["Q100_14dB", "Q100_08dB"]
METHODS = ["Standard", "SC", "SC0", "FullCalRefit", "QScan_Standard"]

reps, infos = [], []
for m in MODELS:
    T = np.asarray(loadmat(ROOT / "cache/inputs/fd2d_20261025/truth" / f"{m}_truth.mat")["target"], float)
    nT = np.linalg.norm(T)
    for s in SCENES:
        for r in range(7):
            L = loadmat(SWEEP / f"{m}_{s}" / f"r{r:03d}.mat", simplify_cells=True)
            best = int(L["info"]["best"])
            out = {"Standard": L["STD"][:, :, 0], "FullCalRefit": L["FCR"][:, :, 0], "QScan_Standard": L["STD"][:, :, best]}
            for lab in ("SC", "SC0"):
                f = ROOT / "cache/solved" / f"{m}_{s}_r{r}_{lab}"
                out[lab] = np.load(f.with_suffix(".npz"))["out"]
                info = json.loads(f.with_suffix(".json").read_text())
                infos.append(dict(model=m, scene=s, r=r, label=lab, converged=info["converged"],
                                  kkt=info["kkt_ratio"], iterations=info["iterations"], seconds=info["seconds"]))
            for meth in METHODS:
                e = np.linalg.norm(np.asarray(out[meth], float) - T)
                reps.append(dict(model=m, scene=s, r=r, method=meth, rel=e / nT, nmse=(e / nT) ** 2))

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
        row = dict(model=m, scene=s, carbonate=m in CARBONATE, primary=s in PRIMARY)
        for meth in METHODS:
            row[f"rel_{meth}"] = val(m, s, meth)
            row[f"nmse_{meth}"] = val(m, s, meth, "nmse")
            row[f"rel0_{meth}"] = val(m, s, meth, rs=[0])
        row["ratio_Standard_SC"] = row["rel_Standard"] / row["rel_SC"]
        row["sc_converged_r1_6"] = all(i["converged"] for i in infos if i["model"] == m and i["scene"] == s
                                       and i["label"] == "SC" and i["r"] >= 1)
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


prim = [x for x in scenes if x["primary"]]
std_better, sc_better = gate(prim, "Standard", "SC"), gate(prim, "SC", "Standard")
verdict = ("Standard stably better than SC" if std_better["stable"] else
           "SC stably better than Standard" if sc_better["stable"] else "no stable advantage")
conv = [x for x in prim if x["sc_converged_r1_6"]]
desc = {}
for name, rows in {"14dB": [x for x in prim if x["scene"] == "Q100_14dB"],
                   "8dB": [x for x in prim if x["scene"] == "Q100_08dB"],
                   "carbonate": [x for x in prim if x["carbonate"]],
                   "no_carbonate": [x for x in prim if not x["carbonate"]],
                   "Q070_14dB": [x for x in scenes if x["scene"] == "Q070_14dB"]}.items():
    desc[name] = dict(Standard_vs_SC=gate(rows, "Standard", "SC"))
desc["primary_other_methods_mean_rel"] = {k: float(np.mean([x[f"rel_{k}"] for x in prim])) for k in METHODS}
desc["primary_mean_nmse"] = {k: float(np.mean([x[f"nmse_{k}"] for x in prim])) for k in METHODS}
desc["noise_free_r0_mean_rel"] = {k: float(np.mean([x[f"rel0_{k}"] for x in prim])) for k in METHODS}
desc["SC_vs_SC0_primary"] = gate(prim, "SC", "SC0")
cmn.writejson(an / "decisions.json", dict(
    verdict=verdict, standard_better=std_better, sc_better=sc_better,
    sc_solves=len(infos), sc_unconverged=[i for i in infos if not i["converged"]],
    sensitivity_converged_only=dict(n=len(conv), standard_better=gate(conv, "Standard", "SC") if conv else None,
                                    sc_better=gate(conv, "SC", "Standard") if conv else None),
    descriptive=desc))
for x in prim:
    print(f"{x['model']} {x['scene']}: Standard {x['rel_Standard']:.4f} SC {x['rel_SC']:.4f} ratio {x['ratio_Standard_SC']:.3f} "
          f"SC0 {x['rel_SC0']:.4f} FCR {x['rel_FullCalRefit']:.4f} QScan+Std {x['rel_QScan_Standard']:.4f}")
print(json.dumps(dict(verdict=verdict, standard_better=std_better, sc_better=sc_better), indent=1))
