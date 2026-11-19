"""Data for paper figure 8: one gate-A2 gather (M24, 60 Hz, events_aligned, seed 0, 24 dB), c = 0.5 and 2.0. Per-event AVO
amplitude vs sin^2 for the truth and the oracle / uncalibrated / calibrated saturated damping (same c_hat as gate A2).
Amplitude at offset x = r_true(x) <bp O_x, bp T_x> / <bp T_x, bp T_x> over the event window (as step2b.metrics, band 10 Hz - fp)."""
import sys, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / "legacy2b"))
import e2e
from e2e import s

M, fp, scene, seed, snr = "M24", 60.0, "events_aligned", 0, 24; Qm = np.array(e2e.MODELS[M], float); mnum = 24
s.g.RIC = s.spec(fp); s.th.S2 = s.g.RIC ** 2; e2e.SNRS = (snr,)
A2 = json.loads((e2e.ROOT / "results" / "a2_M24_f60.json").read_text()); out = []
for c in (0.5, 2.0):
    ce = [r["c_hat"] for r in A2 if r["method"] == "cal" and r["c"] == c and r["snr"] == snr and r["seed"] == seed and r["scene"] == scene][0]
    D, Cl, Tr, sig, meta, keep = e2e.columns(c * Qm, mnum, scene, fp, [seed])
    res = {"oracle": e2e.solve(D, Cl, sig, keep, c * Qm, "sat")[0], "uncalibrated": e2e.solve(D, Cl, sig, keep, Qm, "sat")[0], "calibrated": e2e.solve(D, Cl, sig, keep, ce * Qm, "sat")[0]}
    bT = s.bp(Tr[:, :, 0], (10.0, fp)); ei = np.round(s.C.event_t / s.DT).astype(int); win = 3.0 / fp
    for j, i0 in enumerate(ei):
        sin_i = np.array([s.g.GEO[i]["p"][i0] * s.C.vp[j] for i in range(s.NX)]); rt = s.g.rpp(j, sin_i); w = np.abs(s.T - s.T[i0]) <= win
        row = dict(c=c, c_hat=ce, event=j, sin2=s.g.SIN2[i0].tolist(), truth=rt.tolist())
        for name, O in res.items():
            bO = s.bp(O[:, :, 0], (10.0, fp)); row[name] = (rt * np.sum(bO[w] * bT[w], 0) / np.sum(bT[w] ** 2, 0)).tolist()
        out.append(row)
    print(c, ce, flush=True)
(e2e.ROOT / "results" / "fig_example.json").write_text(json.dumps(out, indent=1))
