"""Second grid extension of the experiment-12 development (PLAN_ZH.md section 5): SC-a improved monotonically as the
sparsity weight a fell (best a = 1/32 at the edge), so a goes further down and every a in A_K also gets the
decoupled variant (amplitude on c_SR x k, phase fix to c_SR), named SCaK_a<a>_<k>+ph.

usage: dev12_ext2.py K N    worker K of N over (dataset, solve) jobs; writes analysis/dev_ext2/<id>__<name>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib12 as L  # noqa: E402
from lib12 import np, sc11, semi, lib  # noqa: E402

A_ONE = [0.015625, 0.0078125, 0.00390625]
A_K, KS = [0.03125, 0.015625, 0.0078125, 0.00390625], [1.3, 1.6]
OUT = L.ROOT / "analysis" / "dev_ext2"


def jobs():
    ids = sorted(p.stem for p in (L.ROOT / "analysis" / "dev").glob("H_*.json"))
    js = [(d, a, 1.0) for d in ids for a in A_ONE] + [(d, a, k) for d in ids for a in A_K for k in KS]
    return js


def name(a, k):
    return f"SCa_a{a}" if k == 1.0 else f"SCaK_a{a}_{k}+ph"


def main(kw, n):
    OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); prep = {}
    for i, (did, a, k) in enumerate(jobs()):
        path = OUT / f"{did}__{name(a, k)}.json"
        if i % n != kw or path.exists():
            continue
        t0 = time.monotonic(); z = np.load(L.DEV_E10 / f"{did}.npz"); D = z["Input"].astype(float)
        ideal = z["ideal"].astype(float)
        c_sr = json.loads((L.ROOT / "analysis" / "dev" / f"{did}.json").read_text(encoding="utf-8"))["c_sr"]
        if did not in prep:
            prep = {did: sc11.prepare_sc(D, cm)}
        p = prep[did]
        o = lib.solve(cm, Q * c_sr * k, p["wv"], D, p["sigma"], p["K"], p["live"], a, lib.SEL_SC[1])[0]
        if k != 1.0:
            o = sc11.phase_fix(o, cm, Q, c_sr * k, c_sr)
        rec = dict(id=did, name=name(a, k), a=a, k=k, score=L.score(o, ideal), seconds=time.monotonic() - t0)
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, rec["name"], f"{rec['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
