"""Experiment 13 runs (PLAN_ZH.md).

usage: run13.py lines        well paths -> inlines to process (geometry only); writes analysis/lines.json
       run13.py qscan K N    QScan factor on QSCAN_LINES; results/qscan/IL_<L>.json
       run13.py a K N        SR-based methods on every inline; cache/out/IL_<L>_A.npz, results/lines/IL_<L>_A.json
       run13.py b K N        the user's pipeline (needs every QScan line); cache/out/IL_<L>_B.npz, ..._B.json
Workers K of N skip finished work (resumable). Logs print no scores.
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib13 as L  # noqa: E402
from lib13 import np, semi, sc11, lib, tie, wq  # noqa: E402

OUT = L.ROOT / "results"; CACHE = L.ROOT / "cache" / "out"; AN = L.ROOT / "analysis"


def well_logs():
    P = wq.picks(); logs = {n: tie.tie_log(n, P) for n in tie.TIE_WELLS}
    return {n: g for n, g in logs.items() if g["t1"] - g["t0"] >= tie.MIN_SPAN}


def lines():
    logs = well_logs(); h, x, y, il, xl = L.cube(); per = {}
    for n, g in logs.items():
        w, idx = tie.path(g, x, y, il, xl); u = np.unique(idx); per[n] = sorted(set(int(v) for v in il[u]))
        assert np.all((xl[u] >= L.XL0) & (xl[u] <= L.XL1) & ((xl[u] - L.XL0) % L.XSTEP == 0)), n
    allL = sorted(set(sum(per.values(), [])) | {semi.INLINE})
    AN.mkdir(exist_ok=True)
    (AN / "lines.json").write_text(json.dumps(dict(lines=allL, per_well=per), indent=1), encoding="utf-8")
    print(len(allL), "inlines")


def todo():
    return json.loads((AN / "lines.json").read_text(encoding="utf-8"))["lines"]


def qscan(k, n):
    (OUT / "qscan").mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    for i, ln in enumerate(L.QSCAN_LINES):
        path = OUT / "qscan" / f"IL_{ln}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); D = L.read_inline(ln); p = sc11.prepare_sc(D, cm); fq, info = sc11.qscan(D, cm, Q, p)
        path.write_text(json.dumps(dict(inline=ln, c=fq, at_boundary=bool(info["at_boundary"]), grid=info["grid"],
                                        objective=info["objective"], seconds=time.monotonic() - t0), indent=1),
                        encoding="utf-8")
        print(ln, f"{time.monotonic() - t0:.0f}s", flush=True)


def part(tag, k, n):
    (OUT / "lines").mkdir(parents=True, exist_ok=True); CACHE.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); qs = None
    if tag == "B":
        qs = {ln: json.loads((OUT / "qscan" / f"IL_{ln}.json").read_text(encoding="utf-8"))["c"] for ln in L.QSCAN_LINES}
    for i, ln in enumerate(todo()):
        path = OUT / "lines" / f"IL_{ln}_{tag}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); D = L.read_inline(ln); p = sc11.prepare_sc(D, cm); c = L.c_sr(ln)
        if tag == "A":
            outs = L.methods_a(D, cm, Q, c, p); info = dict(c_sr=c)
        else:
            fq = L.c_qscan(ln, qs); outs = L.methods_b(D, cm, Q, fq, p); info = dict(c_qscan=fq)
        rec = dict(inline=ln, **info, live=int(p["live"].sum()), bands={m: L.bands(o) for m, o in outs.items()},
                   seconds=time.monotonic() - t0)
        np.savez_compressed(CACHE / f"IL_{ln}_{tag}.npz", **{m: o.astype(np.float32) for m, o in outs.items()})
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(ln, tag, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "lines":
        lines()
    elif m == "qscan":
        qscan(int(sys.argv[2]), int(sys.argv[3]))
    else:
        part(m.upper(), int(sys.argv[2]), int(sys.argv[3]))
