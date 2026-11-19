"""Descriptive: first-order formula dGA ~ -k (u - 1) (uq_avo_20261118/scripts/theory_k.py) against gate A2's NOISE-FREE signed
per-event G/A ('gac0..2', same operator and lambda on clean data) of the uncalibrated saturated damping, events_aligned.
B/A from the generator's exact Zoeppritz coefficients (legacy2b.gate1.rpp), band weight as step2b.metrics."""
import sys, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE / "legacy2b"))
import gate1 as g  # legacy
import e2e

F = np.fft.rfftfreq(4096, g.DT); B = -np.pi * F
out = []
R = [r for m in ("M24", "M25", "M26", "M27") for fp in (60, 40) for r in json.loads((e2e.ROOT / "results" / ("a2_%s_f%d.json" % (m, fp))).read_text())]
cells = {}
for r in R:
    if r["method"] == "uncal" and r["scene"] == "events_aligned" and r["c"] != 1.0:
        cells.setdefault((r["model"], r["fp"], r["c"], r["snr"]), []).append(r)
for (M, fp, c, snr), rows in cells.items():
    Qm = np.array(e2e.MODELS[M], float); a0, ax = g.a0_ax(Qm); ei = np.round(g.C.event_t / g.DT).astype(int); u = 1 / c
    lam = float(np.exp(np.mean(np.log([r["lam"] for r in rows]))))
    S = F ** 2 / fp ** 2 * np.exp(-F ** 2 / fp ** 2); W = S ** 2 * (np.clip((F - 5) / 5, 0, 1) * np.clip((fp + 8 - F) / 8, 0, 1)) ** 2
    for j, i0 in enumerate(ei):
        s2 = g.SIN2[i0]; axj = ax[i0]; G0 = np.exp(-2 * np.pi * F * a0[i0]); rho = G0 / (G0 + lam); feff = np.sum(W * rho * F) / np.sum(W * rho)
        sin_i = np.array([g.GEO[i]["p"][i0] * g.C.vp[j] for i in range(g.NX)]); co = g.rpp(j, sin_i)
        kappa = np.polyfit(s2, axj, 1)[0]; beta, alpha = np.polyfit(s2, s2 ** 2, 1); b1, b0 = np.polyfit(s2, co, 1); ba = b1 / b0
        k = np.pi * feff * kappa * (1 + ba * beta - ba ** 2 * alpha)
        meas = float(np.mean([r["gac%d" % j] for r in rows]))
        out.append(dict(model=M, fp=fp, c=c, snr=snr, ev=j, BA=float(ba), pred_lin=float(-k * (u - 1)), pred_log=float(-k * np.log(u)), meas=meas))
x = np.array([r["meas"] for r in out])
for key in ("pred_lin", "pred_log"):
    y = np.array([r[key] for r in out])
    print("%s: %d points, measured vs predicted slope %.3f R2 %.3f" % (key, len(x), np.polyfit(y, x, 1)[0], np.corrcoef(x, y)[0, 1] ** 2))
    for j in range(3):
        m = np.array([r["ev"] == j for r in out]); print("   event %d (B/A %.2f): slope %.3f R2 %.3f" % (j, out[j]["BA"], np.polyfit(y[m], x[m], 1)[0], np.corrcoef(x[m], y[m])[0, 1] ** 2))
    for c in (0.5, 0.7, 1.4, 2.0):
        m = np.array([r["c"] == c for r in out]); print("   c %.1f: pred median %+.3f, measured median %+.3f" % (c, np.median(y[m]), np.median(x[m])))
(e2e.ROOT / "results" / "check_k_a2.json").write_text(json.dumps(out, indent=1))
