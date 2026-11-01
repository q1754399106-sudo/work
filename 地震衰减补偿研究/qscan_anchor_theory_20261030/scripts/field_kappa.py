"""Field kappa (PLAN_FIELD_ZH.md): measure the anchor coefficient on real Volve inlines, no truth needed.

usage: field_kappa.py [--smoke]    worker; jobs claimed with results9/<job>.claim, resumable.
Lines: the ten QScan inlines of experiment 13 (lib13.QSCAN_LINES), read as in experiment 13 (lib13.read_inline,
xlines 2100-2500 every 2, 1.2-3.0 s), curtain model Q (semi.curtain_cm), lambda 1/2, sc11.prepare_sc.
Scans (as qprofile.scan: from small v upwards, stop 3 grid points past the running minimum, vertex of the
parabola in v = 1/c through the grid minimum and its neighbours):
  default wavelet window (lib.wavelets, 1.2-1.6 s / 1.4 s): c_a = 0.7, 1.0, 1.4
  deep wavelet window (anchor2.wavelets_at 2.4-2.8 s / 2.6 s): c_a = 0.7, 1.4
-> results9/IL_<inline>__<window>_a<c_a>.json
"""
import os, sys, json, time
from env15 import np, lib, semi, sc11, bm, ROOT
from anchor2 import wavelets_at
from qprofile import scan

sys.path.insert(0, str(ROOT.parent / "field_baselines_volve_20261029" / "scripts"))
import lib13 as L13  # noqa: E402

OUT = ROOT / "results9"; bm.MAXITER = 20000
SCANS = [("w1.4", 0.7), ("w1.4", 1.0), ("w1.4", 1.4), ("w2.6", 0.7), ("w2.6", 1.4)]
WIN = {"w1.4": None, "w2.6": (2.4, 2.8, 2.6)}


def jobs():
    return [(ln, w, a) for ln in L13.QSCAN_LINES for w, a in SCANS]


def run(ln, w, ca, cm, Q):
    D = L13.read_inline(ln); pp = sc11.prepare_sc(D, cm); acm = dict(cm, Q=Q * ca)
    wv = lib.wavelets(D, semi.TW, acm) if WIN[w] is None else wavelets_at(D, semi.TW, acm, *WIN[w])
    t0 = time.monotonic(); res = scan(D, cm, Q, pp, lambda cc: wv)
    return dict(inline=ln, window=w, c_anchor=ca, v_anchor=1 / ca, lam=0.5, live_traces=int(pp["live"].sum()),
                seconds=time.monotonic() - t0, **res)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    if "--smoke" in sys.argv:          # timing and reproduction of experiment 13 on one line (a1.0 only)
        r = run(10141, "w1.4", 1.0, cm, Q); r.pop("J"); r.pop("converged"); print(json.dumps(r, indent=1)); sys.exit()
    for ln, w, ca in jobs():
        job = f"IL_{ln}__{w}_a{ca}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(ln, w, ca, cm, Q)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f} grid {rec['v_grid']:.3f}", "edge" if rec["edge"] else "",
              f"{rec['seconds']:.0f}s", flush=True)
