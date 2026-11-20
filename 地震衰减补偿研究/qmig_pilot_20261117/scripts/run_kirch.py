"""Kirchhoff pilot: gain laws x Q multipliers x noise (exploratory, synthetic). Results in results/kirch.json.

Noise: complex Gaussian, white in the 3-42.5 Hz band, level relative to the mean power of the attenuated data.
LP cutoff and STAB gain limit are tuned per noise level at m = 1 on the target-zone E of realisation 0
(favourable to the two baselines); WIEN and UQ have no tuned parameter (UQ: sig = 0.447 as in UQ-W; UQ25: 0.25).
"""
import sys
from kirch import *

LPS = (14, 18, 22, 26, 30, 36); GLS = (6, 10, 15, 20, 30, 40)
NOISE = (("none", None, 1), ("30dB", 30.0, 2), ("14dB", 14.0, 2))

if __name__ == "__main__":
    t0 = time.time(); d = np.load(CACHE / "data.npz"); Dv, Dr = d["visco"], d["ref"]; T, A = tables()
    iz = np.arange(len(IXR)); tq = float(np.median((T[iz, iz, 40:] / np.maximum(A[iz, iz, 40:], 1e-9))))
    amax = 2 * float(A.max())
    xk0 = np.linspace(0, amax, KBIN)
    R = migrate(Dr, law("NONE", xk0), 1.0, xk0, T, A)
    res = dict(tq=tq, amax=amax, ref_rms=float(R.std())); keep = dict(REF=R)
    print("tables %.0f s, tq = %.1f (effective Q), max path t* = %.4f s" % (time.time() - t0, tq, amax), flush=True)
    psig = (abs(Dv) ** 2).mean()
    for tag, snr, nreal in NOISE:
        res[tag] = {}
        for real in range(nreal):
            rng = np.random.default_rng([SEED, int(snr or 0), real]); D = Dv; Nz = None
            if snr is not None:
                sd = np.sqrt(psig / 10 ** (snr / 10) / 2); g = lambda: sd * (rng.standard_normal(Dv.shape) + 1j * rng.standard_normal(Dv.shape))
                D = Dv + g(); Nz = g(); inv0 = inv_snr_windows(D, Nz)
                In = migrate(Nz, law("NONE", xk0), 1.0, xk0, T, A); Id = migrate(D, law("NONE", xk0), 1.0, xk0, T, A)
                sx, sz = ZONES["all"]; pn = (In[sx, sz] ** 2).mean(); snr_img = max((Id[sx, sz] ** 2).mean() - pn, 1e-3 * pn) / pn
                pd = (abs(D) ** 2).mean(); pz = (abs(Nz) ** 2).mean(); fold = snr_img / (max(pd - pz, 1e-3 * pz) / pz)
                res[tag].setdefault("fold", []).append(float(fold)); inv0 = inv0 / fold
                print(tag, real, "trace snr %.1f dB, image snr %.1f dB, fold %.0f" % (10 * np.log10(snr_img / fold), 10 * np.log10(snr_img), fold), flush=True)
            def run(name, m, par=None):
                xk = np.linspace(0, amax / m, KBIN)
                inv = inv0 if (name in ("WIEN", "UQ") and Nz is not None) else None
                return migrate(D, law(name, xk, par, inv), m, xk, T, A)
            if real == 0:
                tune = {}
                for name, grid in (("LP", LPS), ("STAB", GLS)):
                    sc = [score(run(name, 1.0, p), R)["target"]["E"] for p in grid]
                    tune[name] = dict(grid=list(grid), E=sc, best=grid[int(np.argmax(sc))])
                    print(tag, name, "tuning E(target):", " ".join("%g:%.3f" % z for z in zip(grid, sc)), flush=True)
                res[tag]["tune"] = tune
            out = {}
            for m in MS:
                for name, par in (("NONE", None), ("PH", None), ("LP", tune["LP"]["best"]), ("STAB", tune["STAB"]["best"]),
                                  ("WIEN", None), ("UQ", None), ("UQ25", 0.25)):
                    I = run("UQ" if name == "UQ25" else name, m, par); out["%s_m%g" % (name, m)] = score(I, R)
                    if real == 0:
                        keep["%s_%s_m%g" % (tag, name, m)] = I.astype(np.float32)
                print(tag, real, "m=%g" % m, " ".join("%s %.3f/%.3f" % (n, out["%s_m%g" % (n, m)]["target"]["corr"],
                      out["%s_m%g" % (n, m)]["target"]["E"]) for n in ("NONE", "PH", "LP", "STAB", "WIEN", "UQ", "UQ25")),
                      "[%.0f s]" % (time.time() - t0), flush=True)
            res[tag]["real%d" % real] = out
    json.dump(res, open(ROOT / "results" / "kirch.json", "w"), indent=1)
    np.savez_compressed(CACHE / "kirch_images.npz", **keep)
