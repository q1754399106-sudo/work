"""Scoring of PLAN_ZH.md: Yan 2023 selections for phi = 0, 0.1, ..., 1, offsets against the
archived error-optimal grid point, the Y1 judgment and the registered descriptive quantities.
Reads cache/yan/*.mat (yan_stats.m) and the archived summaries/curves in cache/inputs.
Writes analysis/selections.csv, analysis/offsets.csv and analysis/decisions.json.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
import csv, json  # noqa: E402,E401
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.io import loadmat  # noqa: E402

GRID = np.exp(np.linspace(np.log(.6), np.log(2.2), 13))
PHIS = [round(0.1 * i, 1) for i in range(11)]
MODELS = [f"G{m}" for m in range(1, 9)] + [f"H{m}" for m in range(1, 7)]
SCENES = ["Q100_14dB", "Q100_08dB", "Q070_14dB"]
IN = ROOT / "cache" / "inputs"
summ = pd.concat([pd.read_csv(IN / "qscan_offset_20261025/analysis/offset_summary.csv"),
                  pd.read_csv(IN / "fd2d_20261025/analysis/fd_summary.csv")], ignore_index=True)
curves = pd.concat([pd.read_csv(IN / "qscan_offset_20261025/analysis/offset_curves.csv"),
                    pd.read_csv(IN / "fd2d_20261025/analysis/fd_curves.csv")], ignore_index=True)
gidx = lambda v: int(np.argmin(np.abs(GRID - v))) + 1      # 1-based grid index


def minmax(v):
    lo, hi = v.min(), v.max()
    return np.zeros_like(v) if hi == lo else (v - lo) / (hi - lo)


def pick(E, S, phi, norm):
    if norm == "minmax":
        e, s = minmax(E), minmax(S)
    else:
        e, s = E / E.max(), S / S.max()
    return int(np.argmin((1 - phi) * e + phi * s)) + 1


def mode(xs):                                           # smallest of the most frequent (as MATLAB mode)
    vals, counts = np.unique(xs, return_counts=True)
    return int(vals[counts == counts.max()].min())


sel_rows, trace_rows = [], []
for m in MODELS:
    for s in SCENES:
        for r in range(7):
            Y = loadmat(ROOT / "cache/yan" / f"{m}_{s}_r{r}.mat")
            eps, s12, live = Y["eps_"], Y["s12"], Y["live"].ravel().astype(bool)
            E, S = eps[:, live].sum(1), s12[:, live].sum(1) ** 2
            row = dict(model=m, scene=s, r=r, qscan=int(Y["bestJ"].ravel()[0]))
            for phi in PHIS:
                row[f"yan_{phi}"] = pick(E, S, phi, "minmax")
                row[f"yanmax_{phi}"] = pick(E, S, phi, "max")
            sel_rows.append(row)
            if r >= 1:
                for phi in PHIS:
                    per = [pick(eps[:, j], s12[:, j] ** 2, phi, "minmax") for j in np.flatnonzero(live)]
                    trace_rows.append(dict(model=m, scene=s, r=r, phi=phi, median_index=float(np.median(per)),
                                           p10=float(np.percentile(per, 10)), p90=float(np.percentile(per, 90))))
sel = pd.DataFrame(sel_rows)
an = ROOT / "analysis"
an.mkdir(exist_ok=True)
sel.to_csv(an / "selections.csv", index=False)

off_rows = []
for m in MODELS:
    for s in SCENES:
        sm = summ[(summ.model == m) & (summ.scene == s)].iloc[0]
        # The summaries give q_best_grid as a ratio to Q0, i.e. factor x 0.7 in the Q0x0.7 scene
        # (PLAN_ADDENDUM_ZH.md).
        opt = gidx(sm.q_best_grid / (0.7 if s == "Q070_14dB" else 1.0))
        rs = sel[(sel.model == m) & (sel.scene == s) & (sel.r >= 1)]
        q_mode = mode(rs.qscan.values)
        assert q_mode - opt == int(sm.offset_steps), (m, s, q_mode, opt, sm.offset_steps)
        cv = curves[(curves.model == m) & (curves.scene == s)].set_index("index").mse
        row = dict(model=m, scene=s, generator=m[0], opt_index=opt, qscan_offset=q_mode - opt,
                   qscan_mse=float(cv[q_mode + 1]), opt_mse=float(cv[opt + 1]))
        for phi in PHIS:
            for tag in ("yan", "yanmax"):
                k = mode(rs[f"{tag}_{phi}"].values)
                row[f"{tag}_{phi}_offset"] = k - opt
                row[f"{tag}_{phi}_mse"] = float(cv[k + 1])
        off_rows.append(row)
off = pd.DataFrame(off_rows)
off.to_csv(an / "offsets.csv", index=False)


def stability(ks):
    vals, counts = np.unique(ks, return_counts=True)
    spread, top = int(ks.max() - ks.min()), int(counts.max())
    return dict(offsets=[int(k) for k in ks], spread=spread, modal_count=top, n=len(ks),
                mean_abs=float(np.mean(np.abs(ks))), max_abs=int(np.max(np.abs(ks))),
                both_signs=bool((ks > 0).any() and (ks < 0).any()),
                stable=bool(spread <= 1 and top >= int(np.ceil(11 / 14 * len(ks)))))


p = off[off.scene == "Q100_14dB"]
Y1 = {str(phi): stability(p[f"yan_{phi}_offset"].values) for phi in PHIS}
stable_phis = [phi for phi in PHIS if Y1[str(phi)]["stable"]]
verdict = ("Yan offsets model-dependent for every phi" if not stable_phis
           else f"stable for phi in {stable_phis} (exploratory: chosen on test data)")
desc = dict(qscan_Q0_14dB=stability(p.qscan_offset.values))
for s in SCENES:
    q = off[off.scene == s]
    desc[s] = {"qscan": stability(q.qscan_offset.values),
               **{f"yan_{phi}": stability(q[f"yan_{phi}_offset"].values) for phi in PHIS},
               **{f"yanmax_{phi}": stability(q[f"yanmax_{phi}_offset"].values) for phi in PHIS},
               "mse_mean": {"qscan": float(q.qscan_mse.mean()), "optimum": float(q.opt_mse.mean()),
                            **{f"yan_{phi}": float(q[f"yan_{phi}_mse"].mean()) for phi in PHIS}}}
    for g in ("G", "H"):
        qg = q[q.generator == g]
        desc[s][f"by_generator_{g}"] = {"qscan": stability(qg.qscan_offset.values),
                                        **{f"yan_{phi}": stability(qg[f"yan_{phi}_offset"].values) for phi in PHIS}}
tr = pd.DataFrame(trace_rows)
tr.to_csv(an / "per_trace.csv", index=False)
json.dump(dict(Y1=Y1, Y1_verdict=verdict, stable_phis=stable_phis, descriptive=desc),
          open(an / "decisions.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print("QScan offsets reproduced for all 42 model-scenes")
print("Q0/14 dB QScan:", desc["qscan_Q0_14dB"])
for phi in PHIS:
    y = Y1[str(phi)]
    print(f"phi {phi}: offsets {y['offsets']} spread {y['spread']} modal {y['modal_count']}/14 "
          f"mean|k| {y['mean_abs']:.2f} stable {y['stable']}")
print("Y1:", verdict)
