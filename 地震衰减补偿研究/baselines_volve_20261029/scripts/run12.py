"""Experiment 12 held-out runs (PLAN_ZH.md).

usage: run12.py K N   worker K of N over the held-out datasets (resumable); results/synth/<id>.json
T group only (experiment-8 recipe, F-11 T2 log reflectivity), structure line STRUCT, noise line NOISE[0] in the data.
Baseline parameters come from analysis/dev_selection.json (frozen with the plan).
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib12 as L  # noqa: E402
from lib12 import np, sc11, semi, sq, gen, lib  # noqa: E402

STRUCT, NOISE = 10283, [10107, 10235, 10319]
C_T = [0.65, 0.85, 1.05, 1.25, 1.5]
SEEDS_T = [9701, 9702, 9703, 9704]
NBOOT = 200
OUT = L.ROOT / "results"; CACHE = L.ROOT / "cache"
SEL = json.loads((L.ROOT / "analysis" / "dev_selection.json").read_text(encoding="utf-8"))


def datasets():
    return [(f"T_c{c}_s{s}", c, s) for c in C_T for s in SEEDS_T]


def outputs(D, cm, Q, c_sr, c_true):
    """Methods on one section (PLAN_ZH.md section 3)."""
    p = sc11.prepare_sc(D, cm); wv = p["wv"]
    out = dict(Input=D)
    out["SRSCa_1.3"] = sc11.sca(D, cm, Q, c_sr * 1.3, p)
    out["SRSCa_1.3+ph"] = sc11.phase_fix(out["SRSCa_1.3"], cm, Q, c_sr * 1.3, c_sr)
    out["SRSCa"] = sc11.sca(D, cm, Q, c_sr, p)
    a = SEL["B2"]["a"]
    out["B2_SCa"] = out["SRSCa"] if a == lib.SEL_SC[0] else \
        lib.solve(cm, Q * c_sr, wv, D, p["sigma"], p["K"], p["live"], a, lib.SEL_SC[1])[0]
    if SEL["SCaK"]:
        s = SEL["SCaK"]; kk = s["k"]
        out["SCaK"] = sc11.phase_fix(lib.solve(cm, Q * c_sr * kk, wv, D, p["sigma"], p["K"], p["live"], s["a"],
                                               lib.SEL_SC[1])[0], cm, Q, c_sr * kk, c_sr)
    out["B1_IQ"] = L.iq(D, cm, Q, c_sr, c_sr, SEL["B1"]["g_db"])
    out["B3_VDB"] = L.tvw(L.iq(D, cm, Q, None, c_sr, 0), wv["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"])
    out["IQd"] = L.iq(D, cm, Q, c_sr * SEL["IQd"]["k"], c_sr, SEL["IQd"]["g_db"])
    k = SEL["SRSCa_k"]["k"]
    if k != 1.3:
        out[f"SRSCa_{k}+ph"] = sc11.phase_fix(sc11.sca(D, cm, Q, c_sr * k, p), cm, Q, c_sr * k, c_sr)
    out["SCa_true"] = sc11.sca(D, cm, Q, c_true, p)
    return out


def main(k, n):
    (OUT / "synth").mkdir(parents=True, exist_ok=True); ctx = gen.prepare(STRUCT, NOISE)
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float); f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows()
    T = sq.mw_table(cm, Q, f, wins)
    for i, (did, c, seed) in enumerate(datasets()):
        path = OUT / "synth" / f"{did}.json"
        if i % n != k or path.exists():
            continue
        t0 = time.monotonic(); g = gen.signal(ctx, c, seed)
        N0 = g["level"] * ctx["noise"][0]; D = g["x"] + N0; f_true = semi.true_band(g["x"], N0)
        ideal = semi.model(g["R"], None, g["S"] * g["scale"] * semi.cut(g["f"], f_true), g["f"], float(cm["fref"]))
        sr = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=NBOOT, seed=seed)
        outs = outputs(D, cm, Q, sr["c"], c)
        rec = dict(id=did, c_true=c, seed=seed, level=g["level"], lines=ctx["lines"], f_true=f_true, sr=sr,
                   scores={m: L.score(o, ideal) for m, o in outs.items()}, seconds=time.monotonic() - t0)
        CACHE.joinpath("out").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(CACHE / "out" / f"{did}.npz", ideal=ideal.astype(np.float32),
                            **{m: o.astype(np.float32) for m, o in outs.items()})
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
