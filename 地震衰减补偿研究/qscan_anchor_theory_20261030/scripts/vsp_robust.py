"""Robust VSP interval Q from the archived per-level spectral ratios (PLAN_VSP_ZH.md) -> analysis/vsp_robust.json.

Inputs (archived, not reprocessed): field_vsp_volve_f11_20261026/analysis/posthoc_vsp_q.json (F-11 T2: segments RS,
NI: delta_tau one-way, delta_tstar per level vs the reference level, level_r2, fb) and
field_vsp_volve_20261025/analysis/vsp_q.json + vsp_spectra.mat (F-15A: 190 levels vs level 191, fb).
Estimators of 1/Q = slope of delta_tstar on delta_tau: OLS (reproduces the archive) and Theil-Sen (primary);
90% intervals by a moving-block bootstrap over levels in depth order (block 10, 2000 draws, seed 20261109).
Window-equivalent c (analysis window 1.2-3.0 s TWT, curtain model 78.5 above 2.041 s and 135.18 below):
  F-11 T2: c = dt*_model(1.2-3.0) / [L2 / Q_RS + L3 / Q_NI], RS and NI bootstrapped independently;
  F-15A:   over the part of the window its levels cover (1.2 s to twice its deepest first arrival):
           c = dt*_model(1.2-T) / [(T - 1.2) / Q_int].
"""
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.stats import theilslopes

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent; RNG = np.random.default_rng(20261109); NB, BLOCK = 2000, 10
QM2, QM3, TB = 78.50913615296508, 135.17807648006433, 2.04091


def slope(x, y, how):
    return float(np.polyfit(x, y, 1)[0]) if how == "ols" else float(theilslopes(y, x)[0])


def boot_slopes(x, y, how):
    n = len(x); order = np.argsort(x); x, y = x[order], y[order]; out = []
    for _ in range(NB):
        st = RNG.integers(0, n - BLOCK + 1, int(np.ceil(n / BLOCK))); i = np.concatenate([np.arange(s, s + BLOCK) for s in st])[:n]
        out.append(slope(x[i], y[i], how))
    return np.array(out)


def dtm(t0, t1):
    """Model one-way... two-way differential t* of the curtain model between TWT t0 and t1."""
    up = max(0.0, min(t1, TB) - t0); lo = max(0.0, t1 - max(t0, TB)); return up / QM2 + lo / QM3, up, lo


def summary(x, y):
    o = {}
    for how in ("ols", "theilsen"):
        s = slope(x, y, how); b = boot_slopes(x, y, how)
        o[how] = dict(invQ=s, Q=1 / s if s > 0 else float("inf"), invQ_ci=[float(np.percentile(b, 5)), float(np.percentile(b, 95))], boot=b)
    return o


if __name__ == "__main__":
    f11 = json.loads((BASE / "field_vsp_volve_f11_20261026/analysis/posthoc_vsp_q.json").read_text())
    f15 = json.loads((BASE / "field_vsp_volve_20261025/analysis/vsp_q.json").read_text())
    sp15 = loadmat(BASE / "field_vsp_volve_20261025/analysis/vsp_spectra.mat", simplify_cells=True)
    out = {}
    seg = {s: summary(np.array(f11[s]["delta_tau"], float), np.array(f11[s]["delta_tstar"], float)) for s in ("RS", "NI")}
    d_all, L2, L3 = dtm(1.2, 3.0); out["F11"] = {}
    for how in ("ols", "theilsen"):
        rs, ni = seg["RS"][how], seg["NI"][how]
        c = d_all / (L2 * rs["invQ"] + L3 * ni["invQ"]); cb = d_all / (L2 * rs["boot"] + L3 * ni["boot"])
        cb = cb[np.isfinite(cb) & (cb > 0)]
        out["F11"][how] = dict(Q_RS=rs["Q"], Q_RS_invQ_ci=rs["invQ_ci"], Q_NI=ni["Q"], Q_NI_invQ_ci=ni["invQ_ci"], c_window=float(c),
                              c_window_ci=[float(np.percentile(cb, 5)), float(np.percentile(cb, 95))],
                              c_rs_only=float(d_all / (L2 * rs["invQ"] + L3 / QM3)),
                              c_rs_only_ci=[float(d_all / (L2 * rs["invQ_ci"][1] + L3 / QM3)), float(d_all / (L2 * rs["invQ_ci"][0] + L3 / QM3))])
    fb = np.asarray(sp15["fb"], float); T = float(2 * fb.max()); x = np.array(f15["delta_tau"], float); y = np.array(f15["delta_tstar"], float)
    s15 = summary(x, y); d_cov, _, _ = dtm(1.2, min(T, 3.0)); out["F15A"] = dict(twt_covered=[float(2 * fb.min()), T])
    for how in ("ols", "theilsen"):
        r = s15[how]; seg_len = min(T, 3.0) - 1.2; cb = d_cov / (seg_len * r["boot"]); cb = cb[np.isfinite(cb) & (cb > 0)]
        out["F15A"][how] = dict(Q=r["Q"], invQ_ci=r["invQ_ci"], c_window_covered=float(d_cov / (seg_len * r["invQ"])),
                                c_ci=[float(np.percentile(cb, 5)), float(np.percentile(cb, 95))])
    V = {}
    for w, key in (("F11", "c_window_ci"), ("F15A", "c_ci")):
        lo, hi = out[w]["theilsen"][key]
        V[w] = "excludes 0.65 (supports SR)" if lo > 0.75 else ("supports c~0.65" if hi < 0.95 else "undecided")
    out["verdicts"] = V
    (ROOT / "analysis" / "vsp_robust.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
