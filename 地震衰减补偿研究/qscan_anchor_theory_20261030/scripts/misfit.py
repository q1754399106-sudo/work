"""Misfit ladder (PLAN_MISFIT_ZH.md): semi-synthetic variants that may raise the lambda-independent beta.

usage: misfit.py    worker; jobs claimed with results24/<job>.claim, resumable.
Base recipe as betalam.py (base profile, structure 10227, noise x 0.5, anchor at the truth), seeds 9971/9972 x c 0.7/1.0/1.4.
Variants: R0 (SR only; scans in results23), P1 no dispersion phase in the forward model, P2/P3 source phase -45/+45 deg,
M1 water-layer multiples y(t) = r(t) - 0.25 y(t - 0.108 s) on the reflectivity. Noise level re-calibrated as gen.signal.
Scans at lambda 1/16 and 1/2 (qprofile.scan), SR (sq.estimate). -> results24/<variant>_c<c>_s<seed>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, sq, ROOT
import qprofile as qp

OUT = ROOT / "results24"; bm.MAXITER = 20000
VARIANTS = ("R0", "P1", "P2", "P3", "M1"); SEEDS = (9971, 9972); LAMS = (0.0625, 0.5); TAU, RWB = 0.108, 0.25


def model_amp(R, tstar, S, f):
    """semi.model without the dispersion phase (attenuation exp(-pi f t*) only)."""
    N = 8192; aa = np.clip((120 - f) / 20, 0, 1); aa = 0.5 - 0.5 * np.cos(np.pi * aa); b = -np.pi * f
    off = np.round((semi.TW[:, None] - semi.TF[None, :]) / semi.FINE).astype(int) % N; M = np.empty((len(semi.TW), len(semi.TF)))
    for k in range(len(semi.TF)):
        H = S * aa * np.exp(b * tstar[k]); H[-1] = H[-1].real; M[:, k] = np.fft.irfft(H, N)[off[:, k]]
    return M @ R


def multiples(R):
    lag = int(round(TAU / semi.FINE)); Y = R.copy()
    for i in range(lag, R.shape[0]):
        Y[i] = R[i] - RWB * Y[i - lag]
    return Y


def signal(ctx, c, seed, var):
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE); rng = np.random.default_rng(seed)
    R = semi.reflectivity(ctx["tau"], rng)[0]
    if var == "M1":
        R = multiples(R)
    ph = {"P2": -45.0, "P3": 45.0}.get(var, 0.0); S = semi.source_spectrum(ctx["w_in"], f, ph, cm, c)
    ts = semi.tstar_true(cm, c)
    x = model_amp(R, ts, S, f) if var == "P1" else semi.model(R, ts, S, f, float(cm["fref"]))
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / semi.band_power(x, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + L * N0, semi.TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in gen.LEVEL_GRID]
    return x, float(gen.LEVEL_GRID[int(np.argmin(err))])


def run(var, c, s, base):
    cm = dict(base["cm"]); Q = np.asarray(cm["Q"], float); N0 = base["noise"][0]
    x, level = signal(base, c, s, var); D = x + qp.MULT * level * N0
    rec = dict(variant=var, id=f"c{c}_s{s}", u=1 / c, level=level, c_sr=float(sq.estimate(D, cm, Q)["c"]), scans={})
    if var == "R0":                                   # check: this generator reproduces gen.signal (the betalam data)
        g = gen.signal(base, c, s); D0 = g["x"] + qp.MULT * g["level"] * N0
        rec["check_max_abs_diff_rel"] = float(np.max(np.abs(D - D0)) / np.max(np.abs(D0)))
    else:
        pp = sc11.prepare_sc(D, cm); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * c))
        for lam in LAMS:
            qp.LAM = lam; t0 = time.monotonic(); r = qp.scan(D, cm, Q, pp, lambda cc: wv); r.pop("J"); r.pop("converged")
            rec["scans"][str(lam)] = dict(seconds=time.monotonic() - t0, **r)
    return rec


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); base = None
    for var in VARIANTS:
        for s in SEEDS:
            for c in qp.CS:
                job = f"{var}_c{c}_s{s}"; path = OUT / f"{job}.json"
                if path.exists():
                    continue
                try:
                    os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
                except FileExistsError:
                    continue
                base = base or gen.prepare(qp.STRUCT, qp.NOISE); rec = run(var, c, s, base)
                path.write_text(json.dumps(rec), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
                print(job, "c_sr %.3f" % rec["c_sr"], {k: round(v["v_sel"], 3) for k, v in rec["scans"].items()}, flush=True)
