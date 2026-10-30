"""Experiment 16 runs (PLAN_ZH.md).

usage: run16.py smoke            ZA_c0.75_r0 (development structure A), timing and code check
       run16.py test K N         worker K of N over the 12 sealed datasets (structures D, E, F), resumable
results/<stage>/<id>.json; outputs cache/out/<id>.npz (local).
"""
import sys, json, time
from lib16 import np, lib, semi, sq, sc11, ROOT, TEST, SMOKE, load, methods, score, masks, local, band_last

NBOOT = 0


def run(ids, stage, k=0, n=1):
    out = ROOT / "results" / stage; out.mkdir(parents=True, exist_ok=True); (ROOT / "cache" / "out").mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows()
    T = sq.mw_table(cm, Q, f, wins)
    for i, did in enumerate(ids):
        path = out / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); arr, meta = load(did); D = arr["D"]
        sr = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=NBOOT, seed=meta["rc_seed"])
        p = sc11.prepare_sc(D, cm); c_qs, qi = sc11.qscan(D, cm, Q, p); t1 = time.monotonic(); wv = p["wv"]
        outs = methods(D, cm, Q, float(sr["c"]), float(c_qs), meta["c_eff_cum"], p); mk = masks(arr, arr["gas"])
        rec = dict(id=did, structure=meta["structure"], c0=meta["c0"], r=meta["r"], level=meta["level"], f_true=meta["f_true"],
                   c_eff_cum=meta["c_eff_cum"], c_eff_dif=meta["c_eff_dif"], c_sr=float(sr["c"]), qscan=float(c_qs),
                   qscan_boundary=bool(qi["at_boundary"]), fs=wv["fmax_shallow"], fd=wv["fmax_target"],
                   fs_last=band_last(wv["f"], wv["snr_shallow"], wv["f_peak"]),
                   scores={m: score(o, arr["ideal"]) for m, o in outs.items()},
                   local={m: local(o, arr["ideal"], mk) for m, o in outs.items()},
                   seconds_estimates=t1 - t0, seconds=time.monotonic() - t0)
        np.savez_compressed(ROOT / "cache" / "out" / f"{did}.npz", **{m: o.astype(np.float32) for m, o in outs.items()})
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "smoke":
        run(SMOKE, "smoke")
    else:
        run(TEST, "test", int(sys.argv[2]), int(sys.argv[3]))
