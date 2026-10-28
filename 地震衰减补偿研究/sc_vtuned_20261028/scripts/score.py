"""Scoring (PLAN_ZH.md "判定", "只描述"). Standard and QScan + Standard: archived errors from
vfamily_20261027/analysis/rows.csv (parameters unchanged). SC_tuned / QScan + SC_tuned from cache/eval.
Writes analysis/rows.csv, analysis/decisions.json.
"""
import sys, csv, json, hashlib  # noqa: E401
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
from run import VF, EVAL_SCENES, EVAL_ERRORS, NOISES, tag  # noqa: E402
import numpy as np  # noqa: E402
from scipy.io import loadmat  # noqa: E402

sel = json.loads((ROOT / "analysis/selection.json").read_text(encoding="utf-8"))
assert sel["frozen"] is True
P = (sel["primary"]["a"], sel["primary"]["mu_relative"])
S = (sel["secondary"]["a"], sel["secondary"]["mu_relative"])
hashes = json.loads((VF / "logs/MODEL_HASHES.json").read_text(encoding="utf-8"))
old = {(r["scene"], r["noise"], int(r["q_error"])): r for r in csv.DictReader(open(VF / "analysis/rows.csv", encoding="utf-8"))}

rows = []
for s in EVAL_SCENES:
    tp = VF / "truth" / f"{s}_truth.mat"
    assert hashlib.sha256(tp.read_bytes()).hexdigest() == hashes[f"truth/{s}_truth.mat"], s
    tr = loadmat(tp, simplify_cells=True)["tr"]
    for n in NOISES:
        T = np.asarray(next(t for t in tr if t["name"] == n)["target"], float)
        nT = max(np.linalg.norm(T), 1e-30)
        for e in EVAL_ERRORS:
            o = old[(s, n, e)]
            row = dict(cls=s[0], scene=s, noise=n, q_error=e, qscan_factor=float(o["qscan_factor"]))
            for name, pars in (("", P), ("_secondary", S)):
                if name and pars == P:
                    continue
                for label, key in (("SC", "SC_tuned"), ("QSC", "QScan_SC_tuned")):
                    f = ROOT / "cache/eval" / f"{s}_{n}_Q{e:+d}_{label}_{tag(*pars)}"
                    info = json.loads(f.with_suffix(".json").read_text(encoding="utf-8"))
                    X = np.load(f.with_suffix(".npz"))["out"]
                    row[f"rel_{key}{name}"] = float(np.linalg.norm(X - T) / nT)
                    row[f"converged_{key}{name}"] = bool(info["converged"])
            for k in ("SC", "QScan_SC", "Standard", "QScan_Standard", "QScanRefit"):
                row[f"rel_{k}_archived" if k in ("SC", "QScan_SC") else f"rel_{k}"] = float(o[f"rel_{k}"])
            rows.append(row)

an = ROOT / "analysis"
with open(an / "rows.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)


def gate(rr, a, b):
    """Same three gates as vfamily_20261027/scripts/score.py."""
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


def closure(rr, new, archived, ref):
    r_old = np.mean([x[f"rel_{archived}"] for x in rr]) / np.mean([x[f"rel_{ref}"] for x in rr])
    r_new = np.mean([x[f"rel_{new}"] for x in rr]) / np.mean([x[f"rel_{ref}"] for x in rr])
    return dict(ratio_archived=float(r_old), ratio_tuned=float(r_new),
                gap_closed=float((r_old - r_new) / (r_old - 1)) if r_old != 1 else None)


sel_rows = lambda cls, wrong: [x for x in rows if x["cls"] == cls and (x["q_error"] != 0) == wrong]
conv = lambda rr, key: [x for x in rr if x[f"converged_{key}"]]
dec, concl = {}, {}
for cls in ("W", "X"):
    p1 = both(sel_rows(cls, True), "QScan_SC_tuned", "QScan_Standard")
    p2 = both(sel_rows(cls, False), "SC_tuned", "Standard")
    sides = (p1["stable_side"], p2["stable_side"])
    concl[cls] = ("调参后反超" if "SC side" in sides else
                  "调参后仍输" if sides == ("Standard side", "Standard side") else "差距缩小但未定")
    d = dict(P1_Q_wrong=p1, P2_Q_correct=p2, conclusion=concl[cls],
             P1_converged_only=both(conv(sel_rows(cls, True), "QScan_SC_tuned"), "QScan_SC_tuned", "QScan_Standard"),
             P2_converged_only=both(conv(sel_rows(cls, False), "SC_tuned"), "SC_tuned", "Standard"),
             describe=dict(
                 P1_gap=closure(sel_rows(cls, True), "QScan_SC_tuned", "QScan_SC_archived", "QScan_Standard"),
                 P2_gap=closure(sel_rows(cls, False), "SC_tuned", "SC_archived", "Standard"),
                 tuned_over_archived_SC_Qcorrect=float(np.mean([x["rel_SC_tuned"] / x["rel_SC_archived"] for x in sel_rows(cls, False)])),
                 tuned_over_archived_QSC_Qwrong=float(np.mean([x["rel_QScan_SC_tuned"] / x["rel_QScan_SC_archived"] for x in sel_rows(cls, True)]))))
    if S != P:
        d["secondary_report_only"] = dict(
            P1=both(sel_rows(cls, True), "QScan_SC_tuned_secondary", "QScan_Standard"),
            P2=both(sel_rows(cls, False), "SC_tuned_secondary", "Standard"))
    dec[cls] = d
cs = list(concl.values())
overall = ("V 家族确实有利于标准版的构造" if cs == ["调参后仍输"] * 2 else
           "翻转（至少在该类）来自“谁在这个家族上调过参”" if "调参后反超" in cs else "部分解释")
per_scene = {}
for s in EVAL_SCENES:
    rr = [x for x in rows if x["scene"] == s]
    per_scene[s] = dict(
        QSC_tuned_over_QStd_Qwrong=float(np.mean([x["rel_QScan_SC_tuned"] for x in rr if x["q_error"] != 0]) /
                                         np.mean([x["rel_QScan_Standard"] for x in rr if x["q_error"] != 0])),
        SC_tuned_over_Std_Q0=float(np.mean([x["rel_SC_tuned"] for x in rr if x["q_error"] == 0]) /
                                   np.mean([x["rel_Standard"] for x in rr if x["q_error"] == 0])))
out = dict(selection=dict(primary=sel["primary"], secondary=sel["secondary"]), classes=dec, per_scene=per_scene,
           conclusion=overall)
(an / "decisions.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(dict(conclusion=overall, classes=concl,
                      gates={c: {k: (dec[c][k]["stable_side"], round(dec[c][k]["sc_side"]["mean_ratio"], 3),
                                     dec[c][k]["sc_side"]["wins"], dec[c][k]["sc_side"]["n"],
                                     round(dec[c][k]["sc_side"]["worst_ratio"], 3))
                                 for k in ("P1_Q_wrong", "P2_Q_correct")} for c in dec}), ensure_ascii=False, indent=1))
