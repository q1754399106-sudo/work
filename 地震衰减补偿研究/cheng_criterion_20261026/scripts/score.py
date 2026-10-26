"""Scoring of PLAN_ZH.md: Cheng 2025 shallow/deep consistency (C-A primary; C-B, C-C, C-D) and Ma
2018-type lp (p = 0.1, 0.5, 1) selections, offsets against the archived error-optimal grid point,
the C1 judgment and the registered descriptive quantities. The offset machinery (mode, Q0x0.7
conversion, QScan reproduction assertion, stability) is that of yan_criterion_20261026/scripts/score.py.
Reads cache/cheng/*.mat (cheng_stats.m) and the archived summaries/curves in cache/inputs.
Writes analysis/selections.csv, analysis/offsets.csv and analysis/decisions.json.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
import json  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.io import loadmat  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

GRID = np.exp(np.linspace(np.log(.6), np.log(2.2), 13))
MODELS = [f"G{m}" for m in range(1, 9)] + [f"H{m}" for m in range(1, 7)]
SCENES = ["Q100_14dB", "Q100_08dB", "Q070_14dB"]
CRIT = ["chengA", "chengB", "chengC", "chengD", "lp0.1", "lp0.5", "lp1"]
IN = ROOT / "cache" / "inputs"
summ = pd.concat([pd.read_csv(IN / "qscan_offset_20261025/analysis/offset_summary.csv"),
                  pd.read_csv(IN / "fd2d_20261025/analysis/fd_summary.csv")], ignore_index=True)
curves = pd.concat([pd.read_csv(IN / "qscan_offset_20261025/analysis/offset_curves.csv"),
                    pd.read_csv(IN / "fd2d_20261025/analysis/fd_curves.csv")], ignore_index=True)
gidx = lambda v: int(np.argmin(np.abs(GRID - v))) + 1      # 1-based grid index


def releq(As, Ad):                                      # rows: 13 candidates
    return np.linalg.norm(Ad - As, axis=1) / np.linalg.norm(As, axis=1)


def criteria(Y):
    f = Y["freq"].ravel()
    spec = Y["spec"]                                    # 13 x nfreq x 4 (thirds s/d, halves s/d)
    band, fit = (f >= 5) & (f <= 80), (f >= 10) & (f <= 60)
    s3, d3, s2, d2 = (spec[:, band, w] for w in range(4))
    unit = lambda A: A / np.linalg.norm(A, axis=1, keepdims=True)
    lr = np.log(spec[:, fit, 1] / spec[:, fit, 0])
    slope = np.array([np.polyfit(f[fit], lr[q], 1)[0] for q in range(13)])
    vals = {"chengA": releq(s3, d3), "chengB": releq(unit(s3), unit(d3)), "chengC": releq(s2, d2),
            "chengD": np.abs(slope)}
    P = list(np.round(Y["P"].ravel(), 3))
    for p in (0.1, 0.5, 1.0):
        vals[f"lp{p:g}"] = Y["Lp"][:, P.index(p)]
    # A non-finite value (e.g. an all-zero output) can never be selected.
    return {k: int(np.argmin(np.where(np.isfinite(v), v, np.inf))) + 1 for k, v in vals.items()}


def mode(xs):                                           # smallest of the most frequent (as MATLAB mode)
    vals, counts = np.unique(xs, return_counts=True)
    return int(vals[counts == counts.max()].min())


rows = []
for m in MODELS:
    for s in SCENES:
        for r in range(7):
            Y = loadmat(ROOT / "cache/cheng" / f"{m}_{s}_r{r}.mat")
            rows.append(dict(model=m, scene=s, r=r, qscan=int(Y["bestJ"].ravel()[0]), **criteria(Y)))
sel = pd.DataFrame(rows)
an = ROOT / "analysis"
an.mkdir(exist_ok=True)
sel.to_csv(an / "selections.csv", index=False)

off_rows = []
for m in MODELS:
    for s in SCENES:
        sm = summ[(summ.model == m) & (summ.scene == s)].iloc[0]
        # q_best_grid is a ratio to Q0, i.e. factor x 0.7 in the Q0x0.7 scene (yan_criterion addendum).
        opt = gidx(sm.q_best_grid / (0.7 if s == "Q070_14dB" else 1.0))
        rs = sel[(sel.model == m) & (sel.scene == s) & (sel.r >= 1)]
        q_mode = mode(rs.qscan.values)
        assert q_mode - opt == int(sm.offset_steps), (m, s, q_mode, opt, sm.offset_steps)
        cv = curves[(curves.model == m) & (curves.scene == s)].set_index("index").mse
        row = dict(model=m, scene=s, generator=m[0], opt_index=opt, qscan_index=q_mode,
                   qscan_offset=q_mode - opt, qscan_mse=float(cv[q_mode + 1]), opt_mse=float(cv[opt + 1]))
        for c in CRIT:
            k = mode(rs[c].values)
            row.update({f"{c}_index": k, f"{c}_offset": k - opt, f"{c}_mse": float(cv[k + 1])})
        off_rows.append(row)
off = pd.DataFrame(off_rows)
off.to_csv(an / "offsets.csv", index=False)


def stability(ks, idx, qscan_ks):
    vals, counts = np.unique(ks, return_counts=True)
    spread, top, n = int(ks.max() - ks.min()), int(counts.max()), len(ks)
    rho = spearmanr(ks, qscan_ks).correlation if np.ptp(ks) and np.ptp(qscan_ks) else None
    return dict(offsets=[int(k) for k in ks], spread=spread, modal_count=top, n=n,
                mean_abs=float(np.mean(np.abs(ks))), max_abs=int(np.max(np.abs(ks))),
                both_signs=bool((ks > 0).any() and (ks < 0).any()),
                endpoint_count=int(np.sum((idx == 1) | (idx == 13))),
                spearman_vs_qscan=None if rho is None else float(rho),
                stable=bool(spread <= 1 and top >= int(np.ceil(11 / 14 * n))))


def block(q):
    out = {"qscan": stability(q.qscan_offset.values, q.qscan_index.values, q.qscan_offset.values)}
    for c in CRIT:
        out[c] = stability(q[f"{c}_offset"].values, q[f"{c}_index"].values, q.qscan_offset.values)
    return out


p = off[off.scene == "Q100_14dB"]
C1 = block(p)["chengA"]
if C1["endpoint_count"] >= 7:
    verdict = "degenerate under the registered conventions (modal choice at a grid end in >= 7/14 models)"
elif C1["stable"]:
    verdict = "offset stable across models"
else:
    verdict = "offset model-dependent"
desc = {}
for s in SCENES:
    q = off[off.scene == s]
    desc[s] = dict(all=block(q), by_generator={g: block(q[q.generator == g]) for g in ("G", "H")},
                   mse_mean={"qscan": float(q.qscan_mse.mean()), "optimum": float(q.opt_mse.mean()),
                             **{c: float(q[f"{c}_mse"].mean()) for c in CRIT}})
json.dump(dict(C1=C1, C1_verdict=verdict, descriptive=desc),
          open(an / "decisions.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print("QScan offsets reproduced for all 42 model-scenes")
for c in ["qscan"] + CRIT:
    y = block(p)[c]
    print(f"{c}: offsets {y['offsets']} spread {y['spread']} modal {y['modal_count']}/14 "
          f"mean|k| {y['mean_abs']:.2f} ends {y['endpoint_count']} rho {y['spearman_vs_qscan']}")
print("C1:", verdict)
