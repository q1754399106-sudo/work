"""Experiment 10 runs (PLAN_ZH.md).

usage: run.py synth K N     worker K of N over the held-out semi-synthetic datasets (resumable)
       run.py real          SR on the real inlines, and compensation of inline 10141 with the SR factor
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq, gen  # noqa: E402
from sq import np, lib, bm, semi, TW  # noqa: E402

ROOT = sq.ROOT; OUT = ROOT / "results"; CACHE = ROOT / "cache"
STRUCT_LINE = 10301
NOISE_LINES = [10011, 10031, 10051]
C_TEST = [0.55, 0.7, 0.9, 1.2, 1.4]
SEEDS_H = [9101, 9102, 9103, 9104]
SEEDS_W = [9201, 9202, 9203]
REAL_LINES = list(range(9985, 10370, 16))            # 25 inlines over the cube
NBOOT = 200


def haar_sigma(x):
    d = (x[1::2][: len(x) // 2] - x[0::2][: len(x) // 2]) / np.sqrt(2)
    return float(np.median(np.abs(d)) / 0.6745)


def datasets():
    """(id, group, c, seed, heavy): heavy datasets also run QScan and the compensations."""
    out = [(f"H_c{c}_s{s}", "H", c, s, True) for c in C_TEST for s in SEEDS_H]
    out += [(f"W_c{c}_s{s}", "W", c, s, False) for c in C_TEST for s in SEEDS_W]
    return out


def shift(A, s):
    return np.stack([np.interp(TW + s, TW, a) for a in A.T], axis=1)


def m1(out, ideal):
    """Zero-lag M1 (experiment 8) and M1 allowing one common shift (+-20 ms, 1 ms) in 2.30-2.80 s."""
    w = semi.win(TW, *semi.WIN_LOG); z = float(np.median(semi.colcorr(out[w], ideal[w])))
    lag = max(float(np.median(semi.colcorr(shift(out, s)[w], ideal[w]))) for s in np.arange(-20, 21) / 1000)
    return dict(M1=z, M1_lag=lag)


def compensate(D, cm, Q, factors):
    """IQ30 for each named factor, SC-a with the SR factor, QScan factor; dict of outputs and info."""
    live = np.abs(D).sum(axis=0) > 0; K, _, _ = bm.predictor(D, live)
    sigma = np.array([haar_sigma(D[:, j]) for j in range(D.shape[1])]); wv = lib.wavelets(D, TW, cm)
    out, info = {}, {}
    fq, qi = lib.qscan(cm, Q, wv, D, sigma, K, live); info["qscan"] = dict(factor=fq, **qi)
    factors = dict(factors, QScan=fq)
    for k, c in factors.items():
        out[f"IQ30_{k}"] = semi.inverse_q(D, cm, Q * c, 30)
    out["SCa_SR"], _, _ = lib.solve(cm, Q * factors["SR"], wv, D, sigma, K, live, *lib.SEL_SC)
    return out, info


def synth(k, n):
    todo = [d for i, d in enumerate(datasets()) if i % n == k]
    (OUT / "synth").mkdir(parents=True, exist_ok=True)
    ctx = gen.prepare(STRUCT_LINE, NOISE_LINES); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    assert ctx["lines"] == [STRUCT_LINE] + NOISE_LINES, ctx["lines"]
    f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows(); T = sq.mw_table(cm, Q, f, wins)
    segs = [s for v in __import__("wells").all_segments().values() for s in v]
    for did, grp, c, seed, heavy in todo:
        path = OUT / "synth" / f"{did}.json"
        if path.exists():
            continue
        t0 = time.monotonic(); g = gen.signal(ctx, c, seed, segs if grp == "W" else None)
        rec = dict(id=did, group=grp, c_true=c, seed=seed, level=g["level"], sr={})
        noises = ctx["noise"] if grp == "H" else ctx["noise"][:1]
        for i, N in enumerate(noises):
            rec["sr"][f"n{i}"] = sq.estimate(g["x"] + g["level"] * N, cm, Q, table=T, wins=wins, nboot=NBOOT,
                                             seed=seed * 10 + i)
        if heavy:
            D = g["x"] + g["level"] * ctx["noise"][0]; N0 = g["level"] * ctx["noise"][0]
            f_true = semi.true_band(g["x"], N0)
            ideal = semi.model(g["R"], None, g["S"] * g["scale"] * semi.cut(g["f"], f_true), g["f"], float(cm["fref"]))
            outs, info = compensate(D, cm, Q, dict(model=1.0, SR=rec["sr"]["n0"]["c"], true=c))
            outs["Input"] = D
            rec.update(qscan=info["qscan"], f_true=f_true, scores={m: m1(o, ideal) for m, o in outs.items()})
            (CACHE / "out").mkdir(exist_ok=True)
            np.savez_compressed(CACHE / "out" / f"{did}.npz", ideal=ideal.astype(np.float32),
                                **{m: o.astype(np.float32) for m, o in outs.items()})
        rec["seconds"] = time.monotonic() - t0
        bm.writejson(path, rec)
        print(f"{did} {rec['seconds']:.0f}s", flush=True)        # no scores in the log


def real():
    (OUT / "real").mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows(); T = sq.mw_table(cm, Q, f, wins)
    lines = semi.read_lines(REAL_LINES + [semi.INLINE]); rows = {}
    for L, D in lines.items():
        rows[str(L)] = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=NBOOT, seed=int(L))
        print(L, "done", flush=True)
    bm.writejson(OUT / "real" / "lines.json", rows)
    D = lines[semi.INLINE]; c_sr = rows[str(semi.INLINE)]["c"]
    outs, info = compensate(D, cm, Q, dict(model=1.0, SR=c_sr))
    _, snr_in, _ = lib.coherence_snr(D, TW, 2.0, 3.0)
    spec = {}
    for m, o in dict(outs, Input=D).items():
        fa, snr, P = lib.coherence_snr(o, TW, 2.0, 3.0); spec[m] = dict(snr=snr.tolist(), power=P.tolist())
        spec[m]["band05_95"] = list(lib.spectral(o, o, TW)[1:])
    bm.writejson(OUT / "real" / "inline10141.json", dict(c_sr=c_sr, qscan=info["qscan"], f=fa.tolist(), spectra=spec))
    (CACHE / "out").mkdir(exist_ok=True)
    np.savez_compressed(CACHE / "out" / "real_10141.npz", Input=D.astype(np.float32),
                        **{m: o.astype(np.float32) for m, o in outs.items()})


if __name__ == "__main__":
    if sys.argv[1] == "synth":
        synth(int(sys.argv[2]), int(sys.argv[3]))
    elif sys.argv[1] == "real":
        real()
