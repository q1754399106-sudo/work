"""Analysis-window test (PLAN_TRUNC_ZH.md): does the QScan-W 1/c slope excess come from deep times?

usage: trunc.py    worker; (dataset, scan) jobs claimed with results4/<job>.claim, resumable.
Datasets: experiment-14 recipe (gen.signal on structure line 10227), c 0.7 / 1.0 / 1.4, new seeds 9941-9944,
noise x 0.5. The input wavelet is always lib.wavelets on the full 1.2-3.0 s section (unchanged); only the analysis
window of the inversion and objective is cut to 1.2-T1 s (data rows, dictionary times, predictor K). sigma and live
come from the full section. Scans: QScan-W at lambda 1/16 with T1 = 3.0 / 2.6 / 2.2; fixed-anchor QScan
(default wavelet window, c_a 0.7 / 1.4) at lambda 1/2 with T1 = 2.2.
-> results4/<dataset>__<scan>.json
"""
import os, json, time
from env15 import np, lib, semi, sc11, bm, gen, ROOT
from noise import scan

STRUCT, NOISE = 10227, [10149, 10273, 10355]
CS = [0.7, 1.0, 1.4]; SEEDS = [9941, 9942, 9943, 9944]; MULT = 0.5
T1S_W = [3.0, 2.6, 2.2]; T1_A = 2.2; ANCHORS = [0.7, 1.4]
LAM_W, LAM_A, MAXITER = 2 ** -4, 0.5, 20000
OUT = ROOT / "results4"; bm.MAXITER = MAXITER


def jobs():
    sc = [f"w_t{t}" for t in T1S_W] + [f"a{a}_t{T1_A}" for a in ANCHORS]
    return [(f"c{c}_s{s}", c, s, x) for c in CS for s in SEEDS for x in sc]


def cut(D, cm, pp, t1):
    """Rows 1.2-t1 s of the section, with the model time axis and predictor K cut to match."""
    keep = semi.TW <= t1 + 1e-9; Dt = D[keep]; ct = dict(cm); ct["t"] = np.asarray(cm["t"])[keep]
    K, _, _ = bm.predictor(Dt, pp["live"]); return Dt, ct, dict(pp, K=K)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); ctx = gen.prepare(STRUCT, NOISE); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    for did, c, s, sc in jobs():
        job = f"{did}__{sc}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        t0 = time.monotonic(); g = gen.signal(ctx, c, s); D = g["x"] + MULT * g["level"] * ctx["noise"][0]
        pp = sc11.prepare_sc(D, cm); t1 = float(sc.split("_t")[1]); Dt, ct, pt = cut(D, cm, pp, t1)
        if sc.startswith("w_"):
            lam, ca = LAM_W, None
            sel, diag = scan(Dt, ct, Q, pt, lam, lambda cc: lib.wavelets(D, semi.TW, dict(cm, Q=cc["Q"])))
        else:
            lam, ca = LAM_A, float(sc[1:].split("_")[0]); wv = lib.wavelets(D, semi.TW, dict(cm, Q=Q * ca))
            sel, diag = scan(Dt, ct, Q, pt, lam, lambda cc: wv)
        rec = dict(id=did, scan=sc, noise_mult=MULT, c_true=c, seed=s, level=g["level"], lam=lam, c_anchor=ca, t1=t1,
                   c_sel=sel, diag=diag, seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, sel, f"{rec['seconds']:.0f}s", flush=True)
