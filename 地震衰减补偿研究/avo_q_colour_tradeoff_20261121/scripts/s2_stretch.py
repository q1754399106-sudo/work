"""S2 of PLAN_ZH.md: does NMO stretch cancel the growth of t* with angle? Cached gate-A2 rows only (M24-M27, no new runs).
gac_j = signed clean per-event G/A error. d_j = gac_j(uncal) - gac_j(oracle), paired on (model, fp, c, snr, seed).
S_scene = LS slope of the (snr, seed)-median d_j on ln(1/c); primary c in {0.7, 1.4}; phi = 1 - S_nmo / S_aligned.
usage: s2_stretch.py -> results/s2.json"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; SRC = ROOT.parent / "offset_q_avo_20261118" / "results"
MODELS = ("M24", "M25", "M26", "M27"); FPS = (40, 60); CS = (0.5, 0.7, 1.4, 2.0); PRIM = (0.7, 1.4)


def load():
    rows = []
    for m in MODELS:
        for fp in FPS:
            rows += [r for r in json.loads((SRC / f"a2_{m}_f{fp}.json").read_text()) if r["method"] in ("uncal", "oracle")
                     and r["scene"] in ("events_aligned", "events_nmo")]
    return {(r["model"], int(r["fp"]), r["scene"], r["c"], r["snr"], r["seed"], r["method"]): r for r in rows}


def slope(cs, y):
    x = np.log(1 / np.asarray(cs)); x = x - x.mean(); return float((x * (np.asarray(y) - np.mean(y))).sum() / (x ** 2).sum())


def main():
    R = load(); keys = sorted({k[:6] for k in R}); cells, rows_out = [], []
    d = lambda k, j: R[k + ("uncal",)][f"gac{j}"] - R[k + ("oracle",)][f"gac{j}"]  # noqa: E731
    for m in MODELS:
        for fp in FPS:
            for j in range(3):
                S = {}; Sall = {}; med = {}
                for sc in ("events_aligned", "events_nmo"):
                    med[sc] = {c: float(np.median([d(k, j) for k in keys if k[:4] == (m, fp, sc, c)])) for c in CS}
                    S[sc] = slope(PRIM, [med[sc][c] for c in PRIM]); Sall[sc] = slope(CS, [med[sc][c] for c in CS])
                cells.append(dict(model=m, fp=fp, event=j, S_al=S["events_aligned"], S_nmo=S["events_nmo"], phi=1 - S["events_nmo"] / S["events_aligned"],
                                  phi_all_c=1 - Sall["events_nmo"] / Sall["events_aligned"], med_d=med,
                                  oracle_nmo_gac_med=float(np.median([R[k + ("oracle",)][f"gac{j}"] for k in keys if k[:3] == (m, fp, "events_nmo")]))))
    phic = {(c["model"], c["fp"], c["event"]): c["phi"] for c in cells}
    obs, pred, pred0, ga_obs = [], [], [], []
    for k in keys:
        m, fp, sc, c, snr, seed = k
        if sc != "events_nmo":
            continue
        ka = (m, fp, "events_aligned", c, snr, seed); on = [R[k + ("oracle",)][f"gac{j}"] for j in range(3)]; un = [R[k + ("uncal",)][f"gac{j}"] for j in range(3)]
        da = [d(ka, j) for j in range(3)]; base = np.mean(np.abs(on))
        obs.append(np.mean(np.abs(un)) - base)
        pred.append(np.mean([abs(on[j] + (1 - phic[(m, fp, j)]) * da[j]) for j in range(3)]) - base)
        pred0.append(np.mean([abs(on[j] + da[j]) for j in range(3)]) - base)
        ga_obs.append(R[k + ("uncal",)]["ga"] - R[k + ("oracle",)]["ga"])
        rows_out.append(dict(model=m, fp=fp, c=c, snr=snr, seed=seed, obs=obs[-1], pred=pred[-1], pred0=pred0[-1]))
    phis = np.array([c["phi"] for c in cells]); o, p, p0 = float(np.median(obs)), float(np.median(pred)), float(np.median(pred0))
    out = dict(cells=cells, phi_median=float(np.median(phis)), phi_q=[float(np.quantile(phis, q)) for q in (0.1, 0.25, 0.75, 0.9)],
               phi_all_c_median=float(np.median([c["phi_all_c"] for c in cells])),
               phi_by_event=[float(np.median([c["phi"] for c in cells if c["event"] == j])) for j in range(3)],
               S_al_median=float(np.median([c["S_al"] for c in cells])), S_nmo_median=float(np.median([c["S_nmo"] for c in cells])),
               delta_obs=o, delta_pred=p, delta_pred_phi0=p0, ratio=p / o if o else float("nan"), ratio_phi0=p0 / o if o else float("nan"),
               ga_delta_obs_median=float(np.median(ga_obs)),
               by_c={str(c): dict(obs=float(np.median([r["obs"] for r in rows_out if r["c"] == c])), pred=float(np.median([r["pred"] for r in rows_out if r["c"] == c])),
                                  pred0=float(np.median([r["pred0"] for r in rows_out if r["c"] == c]))) for c in CS})
    out["S2a"] = bool(0.50 <= out["phi_median"] <= 0.85)
    out["S2b"] = bool(np.sign(p) == np.sign(o) and o != 0 and 0.5 <= p / o <= 2.0); out["S2b_phi0_also_fits"] = bool(np.sign(p0) == np.sign(o) and o != 0 and 0.5 <= p0 / o <= 2.0)
    out["S2"] = out["S2a"] and out["S2b"]
    (ROOT / "results").mkdir(exist_ok=True); (ROOT / "results" / "s2.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))


if __name__ == "__main__":
    main()
