"""Experiment 17 runs (PLAN_ZH.md).

usage: run17.py dev K N     worker K of N over the 24 development datasets: SC-a (a = 2^-13) on the mu grid for
                            the SR and the QScan factor, plus IQ70; scores only (results/dev/<id>.json), resumable
       run17.py test K N    worker K of N over the 12 sealed datasets with the frozen choice (CHOICE.json),
                            results/test/<id>.json, outputs cache/out/<id>.npz (local), resumable
"""
import sys, json, time
from lib17 import (np, lib, semi, sq, lib12, ROOT, GRID, DEV, TEST, SEL, A_TUNED, G_IQ, A_ORIG, MU, load, estimates,
                   sc, score)


def setup():
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows()
    return cm, Q, sq.mw_table(cm, Q, f, wins), wins


def dev_one(did, cm, Q, T, wins):
    t0 = time.monotonic(); arr, meta = load(did); D = arr["D"]
    c_sr, c_qs, qi, p = estimates(D, cm, Q, T, wins, meta["rc_seed"])
    S = {"IQ70": score(lib12.iq(D, cm, Q, c_sr, c_sr, G_IQ), arr["ideal"])}
    for route, c in (("SR", c_sr), ("QT", c_qs)):
        for mu in GRID:
            S[f"{route}_mu{mu:g}"] = score(sc(D, cm, Q, c, p, A_TUNED, mu), arr["ideal"])
    return dict(id=did, structure=meta["structure"], c0=meta["c0"], r=meta["r"], c_eff_cum=meta["c_eff_cum"],
                c_sr=c_sr, qscan=c_qs, qscan_boundary=bool(qi["at_boundary"]), scores=S, seconds=time.monotonic() - t0)


def test_one(did, cm, Q, T, wins, choice):
    t0 = time.monotonic(); arr, meta = load(did); D = arr["D"]
    c_sr, c_qs, qi, p = estimates(D, cm, Q, T, wins, meta["rc_seed"])
    m_sr, m_qt = choice["SR"], choice["QT"]
    outs = dict(Input=D, IQ70=lib12.iq(D, cm, Q, c_sr, c_sr, G_IQ),
                VDB=lib12.tvw(lib12.iq(D, cm, Q, None, c_sr, 0), p["wv"]["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"]),
                SR_mu10=sc(D, cm, Q, c_sr, p, A_TUNED, MU), SR_sel=sc(D, cm, Q, c_sr, p, A_TUNED, m_sr),
                QSCa=sc(D, cm, Q, c_qs, p, A_ORIG, MU), QT_mu10=sc(D, cm, Q, c_qs, p, A_TUNED, MU),
                QT_sel=sc(D, cm, Q, c_qs, p, A_TUNED, m_qt))
    rec = dict(id=did, structure=meta["structure"], c0=meta["c0"], r=meta["r"], c_eff_cum=meta["c_eff_cum"], level=meta["level"],
               c_sr=c_sr, qscan=c_qs, qscan_boundary=bool(qi["at_boundary"]), choice=choice,
               fs=p["wv"]["fmax_shallow"], fd=p["wv"]["fmax_target"],
               scores={m: score(o, arr["ideal"]) for m, o in outs.items()}, seconds=time.monotonic() - t0)
    np.savez_compressed(ROOT / "cache" / "out" / f"{did}.npz", **{m: o.astype(np.float32) for m, o in outs.items()})
    return rec


def main(stage, k, n):
    out = ROOT / "results" / stage; out.mkdir(parents=True, exist_ok=True); (ROOT / "cache" / "out").mkdir(parents=True, exist_ok=True)
    ids = DEV if stage == "dev" else TEST; env = setup()
    choice = json.loads((ROOT / "CHOICE.json").read_text(encoding="utf-8")) if stage == "test" else None
    for i, did in enumerate(ids):
        path = out / f"{did}.json"
        if i % n != k or path.exists():
            continue
        rec = dev_one(did, *env) if stage == "dev" else test_one(did, *env, choice)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
