"""Step 1a (no scores): the adjacent-trace coherence SNR of UQ-W against the true SNR of the development sets.

For each development set, noise kind (res / gau) and level (1x, 4x): in the 0.4 s windows centred in 2.3-2.8 s, per band,
  true  = P_signal / P_noise (same taper and smoothing),
  est   = lib.coherence_snr of D (what uqw.snr_rows uses),
  rho_k = cross-power of the noise alone at trace lag k over its power, q_k the same for the signal alone (k = 1..6).
usage: snr_diag.py -> results/snr_diag.json, printed medians over windows and datasets
"""
import json
from env import np, lib, TW, load, DEV, ROOT, BANDS, F5, WIN, noise_of, spectra, cross, starts

rec = []
for i, did in enumerate(DEV):
    arr, meta = load(did)
    for kind in ("res", "gau"):
        for x in (1.0, 4.0):
            N = noise_of(arr, kind, x, 20261114000 + i); D = arr["signal"] + N
            for a in starts(TW):
                if not (2.3 <= a + WIN / 2 <= 2.8):
                    continue
                Fs, Fn = spectra(arr["signal"], TW, a, a + WIN), spectra(N, TW, a, a + WIN)
                Ps, Pn = cross(Fs, 0), cross(Fn, 0); f, est, _ = lib.coherence_snr(D, TW, a, a + WIN)
                for lo, hi in BANDS:
                    b = (F5 >= lo) & (F5 < hi)
                    r = dict(id=did, kind=kind, x=x, t=a + WIN / 2, band=lo, true=float(Ps[b].sum() / Pn[b].sum()), est=float(np.median(est[b])))
                    for k in range(1, 7):
                        r[f"rho{k}"] = float(cross(Fn, k)[b].sum() / Pn[b].sum()); r[f"q{k}"] = float(cross(Fs, k)[b].sum() / Ps[b].sum())
                    rec.append(r)
    print(did, flush=True)
(ROOT / "results" / "snr_diag.json").write_text(json.dumps(rec), encoding="utf-8")
med = lambda rows, key: np.median([r[key] for r in rows])
for kind in ("res", "gau"):
    for x in (1.0, 4.0):
        print(f"noise {kind} x{x:g}   band: true snr | coherence snr | est/true | noise rho1 rho2 rho3 rho4 rho6 | signal q1 q2 q4 q6")
        for lo, hi in BANDS:
            R = [r for r in rec if r["kind"] == kind and r["x"] == x and r["band"] == lo]
            print("   %2d-%2d Hz: %7.3f | %7.3f | x%5.2f | %+.2f %+.2f %+.2f %+.2f %+.2f | %.3f %.3f %.3f %.3f" % (
                lo, hi, med(R, "true"), med(R, "est"), np.median([r["est"] / r["true"] for r in R]),
                med(R, "rho1"), med(R, "rho2"), med(R, "rho3"), med(R, "rho4"), med(R, "rho6"), med(R, "q1"), med(R, "q2"), med(R, "q4"), med(R, "q6")))
