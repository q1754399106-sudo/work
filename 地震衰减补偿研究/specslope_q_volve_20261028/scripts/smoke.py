"""Smoke test before the plan is frozen: one H-type and one W-type dataset on a 41-trace crop, development
seed 7999, c = 1.1 (not a test value). Prints only status and timing, no estimates or scores."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sq, gen, run, wells  # noqa: E402
from sq import np, semi  # noqa: E402

t0 = time.monotonic(); ctx = gen.prepare(run.STRUCT_LINE, run.NOISE_LINES[:1]); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
segs = [s for v in wells.all_segments().values() for s in v]
for segs_ in (None, segs):
    g = gen.signal(ctx, 1.1, 7999, segs_)
    cols = slice(80, 121); semi.CENTRE = 20
    D = (g["x"] + g["level"] * ctx["noise"][0])[:, cols]
    e = sq.estimate(D, cm, Q, nboot=5)
    assert np.isfinite(e["c"]) and len(e["ci90"]) == 2
    if segs_ is None:
        outs, info = run.compensate(D, cm, Q, dict(model=1.0, SR=e["c"], true=1.1))
        assert all(np.all(np.isfinite(o)) for o in outs.values()), "non-finite output"
        ideal = semi.model(g["R"][:, cols], None, g["S"] * g["scale"], g["f"], float(cm["fref"]))
        assert all(np.isfinite(run.m1(o, ideal)["M1"]) for o in outs.values())
    print("ok", "wells" if segs_ else "bootstrap", f"{time.monotonic() - t0:.0f}s", flush=True)
