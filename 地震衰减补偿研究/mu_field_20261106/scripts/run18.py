"""Experiment 18 runs.

usage: run18.py smoke        inline 10141 only: reproduces experiment 13's mu = 10 outputs, times mu = 0.25; no ties, no scores
       run18.py run K N      worker K of N: SC-a with mu = 0.25 on every inline, SR and QScan factors;
                             cache/out/IL_<L>.npz, results/lines/IL_<L>.json. Resumable. Logs print no scores.
"""
import json
import sys
import time

from env18 import L, np, semi, sc11, lib, CACHE13, CACHE, OUT, MU_NEW, lines, qscan_factors, solve


def setup():
    cm = semi.curtain_cm(); return cm, np.asarray(cm["Q"], float), qscan_factors()


def smoke():
    cm, Q, qs = setup(); ln = 10141; D = L.read_inline(ln); t0 = time.monotonic(); p = sc11.prepare_sc(D, cm)
    print("prepare", round(time.monotonic() - t0, 1), "s", flush=True)
    A = np.load(CACHE13 / f"IL_{ln}_A.npz"); B = np.load(CACHE13 / f"IL_{ln}_B.npz")
    print("input identical to experiment 13:", bool(np.array_equal(D.astype(np.float32), A["Input"])), flush=True)
    for name, c, ref in (("SR", L.c_sr(ln), A["B2_SCa"]), ("QT", L.c_qscan(ln, qs), B["QSCa_tuned"])):
        t0 = time.monotonic(); o = solve(D, cm, Q, c, p, 10); t10 = time.monotonic() - t0
        rel = float(np.linalg.norm(o - ref) / np.linalg.norm(ref))
        t0 = time.monotonic(); n = solve(D, cm, Q, c, p, MU_NEW); t25 = time.monotonic() - t0
        print(f"{name}: c {c:.4f}; mu 10 reproduces cache, relative difference {rel:.2e}; seconds mu 10 {t10:.0f}, mu 0.25 {t25:.0f}", flush=True)
        print(f"{name}: output finite {bool(np.all(np.isfinite(n)))}", flush=True)


def run(k, n):
    (OUT / "lines").mkdir(parents=True, exist_ok=True); CACHE.mkdir(parents=True, exist_ok=True)
    cm, Q, qs = setup()
    for i, ln in enumerate(lines()):
        path = OUT / "lines" / f"IL_{ln}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); D = L.read_inline(ln); p = sc11.prepare_sc(D, cm); c, fq = L.c_sr(ln), L.c_qscan(ln, qs)
        outs = dict(SR_mu025=solve(D, cm, Q, c, p, MU_NEW), QT_mu025=solve(D, cm, Q, fq, p, MU_NEW))
        rec = dict(inline=ln, c_sr=c, c_qscan=fq, live=int(p["live"].sum()), bands={m: L.bands(o) for m, o in outs.items()},
                   seconds=time.monotonic() - t0)
        np.savez_compressed(CACHE / f"IL_{ln}.npz", **{m: o.astype(np.float32) for m, o in outs.items()})
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(ln, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "smoke":
        smoke()
    else:
        run(int(sys.argv[2]), int(sys.argv[3]))
