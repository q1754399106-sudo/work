"""Sliding wavelet window on the field (PLAN_SLIDE_FIELD_ZH.md): two-anchor alpha_a as a function of the wavelet-window
centre t0 on the 26 SR inlines.

usage: field_slide.py [--smoke]    worker; jobs claimed with results28/<job>.claim, resumable, several workers allowed.
Lines: field_kappa2.LINES (the 26 inlines of experiment 10, 9985-10369 every 16), read as field_kappa.py
(lib13.read_inline: xlines 2100-2500 every 2, 1.2-3.0 s). Curtain model Q (semi.curtain_cm), lambda 1/2,
sc11.prepare_sc, qprofile-style scan as field_kappa.run (bm.MAXITER 20000).
Wavelet window: [t0 - 0.2, t0 + 0.2] s, absorption taken at t0 (anchor2.wavelets_at), for
t0 = 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6 s; c_a = 0.7 / 1.0 / 1.4. All 546 scans are run fresh with this script,
including the geometries already scanned before (used only for the S0 reproduction check).
--smoke: inline 10141, t0 1.4, c_a 1.0 only; prints v_sel next to results11 (must agree to 0.002).
-> results28/IL_<inline>__t<t0>_a<c_a>.json
"""
import os, sys, json, time
from field_kappa import L13, np, semi, sc11, bm, ROOT
from field_kappa2 import LINES, SRL
from anchor2 import wavelets_at
from qprofile import scan

OUT = ROOT / "results28"; bm.MAXITER = 20000
T0S = (1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6); CAS = (0.7, 1.0, 1.4); HALF = 0.2


def jobs():
    # shallow windows first so that a partial run already covers the curve's lower part
    return [(ln, t0, ca) for t0 in T0S for ln in LINES for ca in CAS]


def run(ln, t0, ca, cm, Q):
    D = L13.read_inline(ln); pp = sc11.prepare_sc(D, cm)
    wv = wavelets_at(D, semi.TW, dict(cm, Q=Q * ca), round(t0 - HALF, 6), round(t0 + HALF, 6), t0)
    t1 = time.monotonic(); res = scan(D, cm, Q, pp, lambda cc: wv)
    return dict(inline=ln, t0=t0, window=[round(t0 - HALF, 6), round(t0 + HALF, 6)], c_anchor=ca, v_anchor=1 / ca,
                lam=0.5, live_traces=int(pp["live"].sum()), c_sr=float(SRL[str(ln)]["c"]),
                fmax_shallow=float(wv["fmax_shallow"]), fmax_target=float(wv["fmax_target"]),
                seconds=time.monotonic() - t1, **res)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    if "--smoke" in sys.argv:
        r = run(10141, 1.4, 1.0, cm, Q); old = json.loads((ROOT / "results11" / "IL_10141__w1.4_a1.0.json").read_text(encoding="utf-8"))
        print(f"smoke v_sel {r['v_sel']:.4f}  results11 {old['v_sel']:.4f}  diff {r['v_sel'] - old['v_sel']:+.4f}  {r['seconds']:.0f}s")
        sys.exit()
    for ln, t0, ca in jobs():
        job = f"IL_{ln}__t{t0:.1f}_a{ca}"; path = OUT / f"{job}.json"
        if path.exists():
            continue
        try:
            os.close(os.open(OUT / f"{job}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            continue
        rec = run(ln, t0, ca, cm, Q)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); (OUT / f"{job}.claim").unlink()
        print(job, f"v {rec['v_sel']:.3f}", "edge" if rec["edge"] else "",
              "" if rec["neighbourhood_converged"] else "unconv", f"{rec['seconds']:.0f}s", flush=True)
