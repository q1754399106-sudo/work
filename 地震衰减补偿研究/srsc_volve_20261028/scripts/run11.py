"""Experiment 11 runs (PLAN_ZH.md).

usage: run11.py synth K N   worker K of N over the held-out datasets (resumable); results/synth/<id>.json
       run11.py real        SR + SC-a and QScan + SC-a on the real inline 10141; results/real/inline10141.json
Groups: T (experiment-8 recipe, F-11 T2 log reflectivity) and L (colour stress: the long Gardner logs of F-12 / F-14
at the centre trace, experiment-9 synth.py), structure line STRUCT, noise line NOISE[0] in the data.
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sc11  # noqa: E402
from sc11 import np, semi, sq, gen, synth, wq, lib, TW  # noqa: E402

STRUCT, NOISE = 10263, [10095, 10215, 10335]
C_T = [0.6, 0.75, 0.95, 1.15, 1.45]
SEEDS_T = [9501, 9502, 9503, 9504]
C_L = [0.75, 1.15]
SEEDS_L = [9601, 9602]
WELLS_L = ["F-12", "F-14"]
K_AMP = [1.3]            # amplitude factor multipliers of the SRSCa_k+ph variants (PLAN_ZH.md)
NBOOT = 200
OUT = sc11.ROOT / "results"; CACHE = sc11.ROOT / "cache"


def datasets():
    out = [(f"T_c{c}_s{s}", "T", None, c, s) for c in C_T for s in SEEDS_T]
    out += [(f"L_{w}_c{c}_s{s}", "L", w, c, s) for w in WELLS_L for c in C_L for s in SEEDS_L]
    return out


def outputs(D, cm, Q, c_sr, c_true=None):
    """Methods on one section (PLAN_ZH.md section 3); returns dict of outputs and info."""
    p = sc11.prepare_sc(D, cm); fq, qi = sc11.qscan(D, cm, Q, p)
    out = dict(Input=D, QSCa=sc11.sca(D, cm, Q, fq, p), SRSCa=sc11.sca(D, cm, Q, c_sr, p),
               IQ30_SR=semi.inverse_q(D, cm, Q * c_sr, 30))
    out["QSCa+ph"] = sc11.phase_fix(out["QSCa"], cm, Q, fq, c_sr)
    for k in K_AMP:
        out[f"SRSCa_{k}+ph"] = sc11.phase_fix(sc11.sca(D, cm, Q, c_sr * k, p), cm, Q, c_sr * k, c_sr)
    if c_true is not None:
        out["SCa_true"] = sc11.sca(D, cm, Q, c_true, p)
    return out, dict(qscan=fq, qscan_boundary=qi["at_boundary"])


def signal_L(ctx, w, c, seed):
    """synth.signal (experiment 9) returning the same fields as gen.signal (x, level, R, scale, S, f)."""
    cm = ctx["cm"]; f = np.fft.rfftfreq(8192, semi.FINE); rng = np.random.default_rng(seed)
    R = synth.reflectivity(ctx, w, rng); S = semi.source_spectrum(ctx["w_in"], f, 0.0, cm, c)
    x = semi.model(R, semi.tstar_true(cm, c), S, f, float(cm["fref"]))
    b = ctx["band"]; k = float(np.sqrt(ctx["Creal"][b].sum() / semi.band_power(x, 1.2, 1.6)[b].sum())); x = x * k
    N0, bq = ctx["noise"][0], ctx["bq"]
    err = [float(np.sqrt(np.mean((np.log(np.maximum(lib.coherence_snr(x + L * N0, TW, 2.0, 3.0)[1][bq], 1e-3))
                                  - np.log(np.maximum(ctx["snr_real"][bq], 1e-3))) ** 2))) for L in gen.LEVEL_GRID]
    return dict(x=x, level=float(gen.LEVEL_GRID[int(np.argmin(err))]), R=R, scale=k, S=S, f=f)


def synth_run(k, n):
    (OUT / "synth").mkdir(parents=True, exist_ok=True); ctx = gen.prepare(STRUCT, NOISE)
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float); f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows()
    T = sq.mw_table(cm, Q, f, wins); P = wq.picks(); W = {w: wq.well(w, P) for w in WELLS_L}
    for i, (did, grp, well, c, seed) in enumerate(datasets()):
        path = OUT / "synth" / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic()
        g = gen.signal(ctx, c, seed) if grp == "T" else signal_L(ctx, W[well], c, seed)
        N0 = g["level"] * ctx["noise"][0]; D = g["x"] + N0; f_true = semi.true_band(g["x"], N0)
        ideal = semi.model(g["R"], None, g["S"] * g["scale"] * semi.cut(g["f"], f_true), g["f"], float(cm["fref"]))
        sr = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=NBOOT, seed=seed)
        outs, info = outputs(D, cm, Q, sr["c"], c)
        rec = dict(id=did, group=grp, well=well, c_true=c, seed=seed, level=g["level"], lines=ctx["lines"], f_true=f_true,
                   sr=sr, **info, scores={m: sc11.m1(o, ideal) for m, o in outs.items()}, seconds=time.monotonic() - t0)
        CACHE.joinpath("out").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(CACHE / "out" / f"{did}.npz", ideal=ideal.astype(np.float32),
                            **{m: o.astype(np.float32) for m, o in outs.items()})
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


def real():
    (OUT / "real").mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    D = semi.read_lines([semi.INLINE])[semi.INLINE]
    L = json.loads((sc11.E10 / "results" / "real" / "lines.json").read_text(encoding="utf-8")); c_sr = L[str(semi.INLINE)]["c"]
    outs, info = outputs(D, cm, Q, c_sr)
    spec = {}
    for m, o in outs.items():
        fa, snr, Pw = lib.coherence_snr(o, TW, 2.0, 3.0); spec[m] = dict(snr=snr.tolist(), power=Pw.tolist(),
                                                                         band05_95=list(lib.spectral(o, o, TW)[1:]))
    bm_out = dict(c_sr=c_sr, **info, f=fa.tolist(), spectra=spec)
    (OUT / "real" / "inline10141.json").write_text(json.dumps(bm_out, indent=1), encoding="utf-8")
    CACHE.mkdir(exist_ok=True)
    np.savez_compressed(CACHE / "real_10141.npz", **{m: o.astype(np.float32) for m, o in outs.items()})


if __name__ == "__main__":
    synth_run(int(sys.argv[2]), int(sys.argv[3])) if sys.argv[1] == "synth" else real()
