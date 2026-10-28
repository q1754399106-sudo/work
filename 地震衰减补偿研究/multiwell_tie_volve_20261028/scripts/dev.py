"""Development of WR on semi-synthetic wells (seen lines of experiment 10; never the real well data).

usage: dev.py K N    worker K of N; writes analysis/dev/<id>.json (resumable)
Grid: wells F-12, F-14; c in C_DEV; seeds SEEDS_DEV; the three noise lines. Per dataset: WR, the same windows without
the log RC (noRC) and SR (experiment 10, 1.2-3.0 s).
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wq, synth  # noqa: E402
from wq import np, semi, gen  # noqa: E402

STRUCT, NOISE = 10301, [10011, 10031, 10051]
WELLS = ["F-12", "F-14"]
C_DEV = [0.6, 0.8, 1.0, 1.25]
SEEDS_DEV = [9301, 9302, 9303]
SM_GRID = [5]
OUT = wq.ROOT / "analysis" / "dev"


def main(k, n):
    OUT.mkdir(parents=True, exist_ok=True)
    ctx = gen.prepare(STRUCT, NOISE); assert ctx["lines"] == [STRUCT] + NOISE, ctx["lines"]
    cm = ctx["cm"]; Q = np.asarray(cm["Q"], float); P = wq.picks()
    jobs = [(w, c, s) for w in WELLS for c in C_DEV for s in SEEDS_DEV]
    W = {w: wq.well(w, P) for w in WELLS}
    for i, (w, c, s) in enumerate(jobs):
        if i % n != k:
            continue
        path = OUT / f"{w}_c{c}_s{s}.json"
        if path.exists():
            continue
        t0 = time.monotonic(); wl = W[w]; wins = wq.well_windows(wl); R = [wq.rc_power(wl, a, b) for a, b in wins]
        g = synth.signal(ctx, wl, c, s); f = np.fft.rfftfreq(wq.sq.NFFT, wq.DT); T = wq.sq.mw_table(cm, Q, f, wins); SW = wq.sq.windows(); TS = wq.sq.mw_table(cm, Q, f, SW)
        rec = dict(well=w, c_true=c, seed=s, level=g["level"], est={})
        for j, N in enumerate(ctx["noise"]):
            secs = synth.sections(g["x"] + g["level"] * N, wins)
            for sm in SM_GRID:
                e = wq.wr(secs, R, wins, cm, Q, sm=sm, nboot=wq.NBOOT, seed=s * 10 + j, table=T)
                rec["est"][f"n{j}_sm{sm}"] = {kk: v for kk, v in e.items() if not kk.startswith("_")}
            e = wq.wr(secs, [np.ones_like(R[0])] * len(wins), wins, cm, Q, nboot=wq.NBOOT, seed=s * 10 + j, table=T)
            rec["est"][f"n{j}_noRC"] = {kk: v for kk, v in e.items() if not kk.startswith("_")}
            e = wq.sq.estimate(g["x"] + g["level"] * N, cm, Q, table=TS, wins=SW)
            rec["est"][f"n{j}_SR"] = {kk: v for kk, v in e.items() if not kk.startswith("_")}
        rec["seconds"] = time.monotonic() - t0
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(path.name, f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
