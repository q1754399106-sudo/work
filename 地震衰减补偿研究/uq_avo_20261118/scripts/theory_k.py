"""First-order AVO sensitivity to the Q level (paper outline item 1; descriptive, no gate).
For an offset-consistent operator with the model Q, the in-band transfer rho(f) is the same on every trace. If the true path
t* is u times the model's, trace x carries an extra exp[-pi f (u - 1) a(x)], so
    ln Gamma(x) ~ const - pi f_eff (u - 1) a_x(x),   f_eff = int W rho f df / int W rho df,   W = S^2 (x band-pass^2)
and with a_x ~ kappa sin^2(theta) per event:   dGA ~ -k (u - 1),   k = pi f_eff kappa x [1 + (B/A) beta - (B/A)^2 alpha],
where alpha, beta are the intercept and slope of sin^4 regressed on sin^2 over the offsets (the Shuey fit of co x Gamma picks up
a B sin^4 term). k_pred0 = the bracket left out.
Checks: (1) k against the closed-form-Gamma slope d dGA / d ln u at u = 1 (sigma-0 offset-consistent operator, ua.op_ua),
many Q models, both wavelets, three lambdas; (2) the prediction |k (u - 1)| against gate A2's measured G/A error of the
UNCALIBRATED saturated damping (offset_q_avo_20261118, noisy, per-event, band-passed, residual-principle lambda).
usage: theory_k.py   -> results/theory_k.json"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate1 as g, ua
from core import ray, C

ROOT = Path(__file__).resolve().parents[1]; A2DIR = ROOT.parent / "offset_q_avo_20261118" / "results"
QS = dict(g.MODELS, M16=[45, 60, 35, 80], M17=[90, 35, 55, 45], M18=[32, 50, 90, 60], M19=[65, 85, 30, 40], M20=[50, 30, 75, 60],
          M21=[38, 95, 45, 70], M22=[80, 55, 40, 30], M23=[60, 40, 110, 50], M24=[42, 70, 55, 95], M25=[75, 45, 30, 65],
          M26=[35, 40, 85, 50], M27=[95, 60, 45, 35])


def bandw(fp):
    b0, b1 = 10.0, fp; return (np.clip((g.F - b0 + 5) / 5, 0, 1) * np.clip((b1 + 8 - g.F) / 8, 0, 1)) ** 2


def k_pred(Q, fp, lam, band=False, full=True):
    a_ev, s2, co = g.paths(np.asarray(Q, float)); a0 = ray(C.event_t, 0.0, np.asarray(Q, float))[0]
    W = g.spec(fp) ** 2 * (bandw(fp) if band else 1.0); out = []
    for j in range(3):
        G0 = np.abs(g.G(a0[j])) ** 2; rho = G0 / (G0 + lam); feff = np.sum(W * rho * g.F) / np.sum(W * rho)
        kappa = np.polyfit(s2[j], a_ev[j] - a0[j], 1)[0]; beta, alpha = np.polyfit(s2[j], s2[j] ** 2, 1)
        b1, b0 = np.polyfit(s2[j], co[j], 1); ba = b1 / b0
        out.append(np.pi * feff * kappa * ((1 + ba * beta - ba ** 2 * alpha) if full else 1.0))
    return np.array(out)


def k_num(Q, fp, lam, e=0.02):
    Q = np.asarray(Q, float); a_ev, s2, co = g.paths(Q); a0 = np.array([ray(C.event_t, 0.0, Q)[0]] * g.NX).T
    H = ua.op_ua(a0, a_ev, 0.0, lam); S2 = g.spec(fp) ** 2
    d = lambda u, j: g.dga(g.gamma_closed(a_ev[j], u, H[j], S2), s2[j], co[j])
    return np.array([-(d(np.exp(e), j) - d(np.exp(-e), j)) / (2 * e) for j in range(3)])


out = {"check1": [], "check2": []}
for M, Q in QS.items():
    for fp in g.FPS:
        for lam in (1e-4, 10 ** -3.5, 1e-3):
            kp, kp0, kn = k_pred(Q, fp, lam), k_pred(Q, fp, lam, full=False), k_num(Q, fp, lam)
            for j in range(3):
                out["check1"].append(dict(model=M, fp=fp, lam=lam, ev=j, k_pred=float(kp[j]), k_pred0=float(kp0[j]), k_num=float(kn[j])))
for evs, tag in (((0, 1, 2), "all"), ((1, 2), "ev12"), ((0,), "ev0")):
    rr = [r for r in out["check1"] if r["ev"] in evs]; x = np.array([r["k_num"] for r in rr])
    for key in ("k_pred", "k_pred0"):
        y = np.array([r[key] for r in rr]); out["check1_%s_%s" % (tag, key)] = dict(slope=float(np.polyfit(x, y, 1)[0]), R2=float(np.corrcoef(x, y)[0, 1] ** 2), median_ratio=float(np.median(y / x)))
        print("check 1 [%s] %s: %d points, k_num %.2f..%.2f, pred vs num slope %.3f R2 %.3f median pred/num %.2f" % (tag, key, len(x), x.min(), x.max(), *out["check1_%s_%s" % (tag, key)].values()))

# check 2: gate A2, uncalibrated saturated damping, aligned scenes, cell = (model, fp, scene, c, snr)
R = [r for m in ("M24", "M25", "M26", "M27") for fp in (60, 40) for r in json.loads((A2DIR / ("a2_%s_f%d.json" % (m, fp))).read_text())]
cells = {}
for r in R:
    if r["method"] in ("uncal", "oracle") and r["scene"] != "events_nmo":
        cells.setdefault((r["model"], r["fp"], r["scene"], r["c"], r["snr"], r["method"]), []).append((r["ga"], r["lam"]))
for (M, fp, sc, c, snr, meth), v in cells.items():
    if meth != "uncal" or c == 1.0:
        continue
    lam = float(np.exp(np.mean(np.log([q[1] for q in v])))); u = 1.0 / c  # truth Q = c x model Q -> path t* = model / c
    kp = k_pred(QS[M], fp, lam, band=True); ga = float(np.mean([q[0] for q in v])); base = float(np.mean([q[0] for q in cells[(M, fp, sc, c, snr, "oracle")]]))
    out["check2"].append(dict(model=M, fp=fp, scene=sc, c=c, snr=snr, lam=lam, ga=ga, ga_oracle=base,
                              pred_lin=float(np.mean(np.abs(kp * (u - 1)))), pred_log=float(np.mean(np.abs(kp * np.log(u))))))
for key in ("pred_lin", "pred_log"):
    p = np.array([r[key] for r in out["check2"]]); m = np.array([r["ga"] for r in out["check2"]])
    out["check2_" + key] = dict(slope=float(np.polyfit(p, m, 1)[0]), R2=float(np.corrcoef(p, m)[0, 1] ** 2), median_ratio=float(np.median(m / p)))
    for c in (0.5, 0.7, 1.4, 2.0):
        sel = [r for r in out["check2"] if r["c"] == c]
        out["check2_%s_c%.1f" % (key, c)] = dict(pred=float(np.median([r[key] for r in sel])), meas=float(np.median([r["ga"] for r in sel])))
    print("check 2 (%s): %d cells, measured vs predicted slope %.3f R2 %.3f median ratio %.2f" % (key, len(p), out["check2_" + key]["slope"], out["check2_" + key]["R2"], out["check2_" + key]["median_ratio"]))
    print("   by c (pred / measured medians):", "  ".join("c %.1f %.3f/%.3f" % (c, out["check2_%s_c%.1f" % (key, c)]["pred"], out["check2_%s_c%.1f" % (key, c)]["meas"]) for c in (0.5, 0.7, 1.4, 2.0)))
(ROOT / "results" / "theory_k.json").write_text(json.dumps(out, indent=1))
