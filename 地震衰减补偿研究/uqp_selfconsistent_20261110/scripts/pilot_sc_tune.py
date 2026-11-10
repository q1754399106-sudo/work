"""Fair SC-a control at 4 x noise (exploratory, NOTE_ZH.md section 11): select (a, mu_relative) of SC-a on 8 of the 24
development datasets (4 x noise, SR route), then score the selected pair on the other 16. UQ-MW has no parameter
tuned on this data, so SC-a gets the advantage of a selection on the same noise level.

  select: Z{A,B,C,D}_c{0.75,1.05}_r0;  evaluate: the other 16 (all r1 and structures E, F)
  grid:   a in {2^-9, 2^-5, 2^-3, 2^-1} x mu_relative in {0.25, 2, 10, 40}  (2^-13 / 0.25 is pilot_semi_sc.py;
          mu 40 added after the first results put the best mu at the old grid edge 10)

usage: pilot_sc_tune.py select K N | best | eval K N A MU
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402,F401
from uqp import np, lib, semi, ROOT  # noqa: E402
from uqw import load  # noqa: E402
from lib17 import DEV, score, sc11  # noqa: E402

OUT = ROOT / "results" / "pilot_sc_tune"
SELECT = [f"Z{k}_c{c}_r0" for k in "ABCD" for c in (0.75, 1.05)]; EVAL = [d for d in DEV if d not in SELECT]
GRID = [(a, mu) for a in (2.0 ** -9, 2.0 ** -5, 2.0 ** -3, 2.0 ** -1) for mu in (0.25, 2.0, 10.0, 40.0)]


def solve(did, a, mu, cache):
    if did not in cache:
        cache.clear(); arr, meta = load(did); D = arr["signal"] + 4.0 * arr["noise"]
        c = json.loads((ROOT / "results" / "pilot_semi" / f"{did}_nd.json").read_text(encoding="utf-8"))["sr"]["c_sr"]
        cache[did] = (arr, D, c, sc11.prepare_sc(D, semi.curtain_cm()))
    arr, D, c, p = cache[did]; cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); t0 = time.monotonic()
    X, info, _ = lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], a, mu)
    return dict(id=did, a=a, mu=mu, c_sr=c, S=float(score(X, arr["ideal"])["S"]), converged=bool(info["converged"]),
                seconds=time.monotonic() - t0)


if __name__ == "__main__":
    mode = sys.argv[1]; OUT.mkdir(parents=True, exist_ok=True); cache = {}
    if mode != "best":
        k, n = int(sys.argv[2]), int(sys.argv[3])
    if mode == "best":                                   # highest mean S over the 8 selection datasets (complete cells only)
        import glob, collections
        G = collections.defaultdict(list)
        for q in glob.glob(str(OUT / "select_*.json")):
            r = json.loads(Path(q).read_text(encoding="utf-8")); G[(r["a"], r["mu"])].append(r["S"])
        full = {key: np.mean(v) for key, v in G.items() if len(v) == len(SELECT)}; b = max(full, key=full.get)
        for key in sorted(full):
            print("a 2^%+.0f mu %-5g S %.3f" % (np.log2(key[0]), key[1], full[key]), file=sys.stderr)
        print(repr(b[0]), repr(b[1])); sys.exit(0)
    if mode == "select":
        jobs = [(d, a, mu) for d in SELECT for a, mu in GRID]
    else:
        a, mu = float(sys.argv[4]), float(sys.argv[5]); jobs = [(d, a, mu) for d in EVAL]
    for j, (d, a, mu) in enumerate(jobs):
        path = OUT / f"{mode}_{d}_a{np.log2(a):+.0f}_mu{mu:g}.json"
        if j % n != k or path.exists():
            continue
        r = solve(d, a, mu, cache); path.write_text(json.dumps(r, indent=1), encoding="utf-8"); print(r, flush=True)
