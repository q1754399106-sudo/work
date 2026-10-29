"""Extra workers for experiment 14, added while run14.py was running (operational only; results are identical).

usage: extra14.py dev|test N_ORIG
Takes the remaining datasets of the stage from the end of the list. A dataset is skipped when its result exists,
when another extra worker has claimed it (results/<stage>/<id>.claim, created exclusively), or when it is among the
next two pending datasets of an original worker (run14.py <stage> k N_ORIG processes i % N_ORIG == k in order), so
the original and extra workers meet in the middle without computing the same dataset twice.
"""
import os, sys, json, time
import run14
from run14 import np, lib, semi, sq, sc11, gen, wq, OUT, LAMS, STRUCT, NOISE, WELLS_L, qscan_w, signal_L

stage, n_orig = sys.argv[1], int(sys.argv[2])
out = OUT / stage
ids = [p.stem for p in run14.dev_sets()] if stage == "dev" else [t[0] for t in run14.test_sets()]


def reserved():
    """The next two pending datasets of each original worker."""
    res = set()
    for k in range(n_orig):
        pend = [d for i, d in enumerate(ids) if i % n_orig == k and not (out / f"{d}.json").exists()]
        res.update(pend[:2])
    return res


def claim(did):
    try:
        os.close(os.open(out / f"{did}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        return True
    except FileExistsError:
        return False


if stage == "dev":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); paths = {p.stem: p for p in run14.dev_sets()}
else:
    ctx = gen.prepare(STRUCT, NOISE); cm = ctx["cm"]; Q = np.asarray(cm["Q"], float)
    f = np.fft.rfftfreq(sq.NFFT, sq.DT); wins = sq.windows(); T = sq.mw_table(cm, Q, f, wins)
    P = wq.picks(); W = {w: wq.well(w, P) for w in WELLS_L}; sets = {t[0]: t for t in run14.test_sets()}

while True:
    todo = [d for d in reversed(ids) if not (out / f"{d}.json").exists() and not (out / f"{d}.claim").exists()
            and d not in reserved()]
    if not todo:
        break
    did = todo[0]
    if not claim(did):
        continue
    t0 = time.monotonic()
    if stage == "dev":
        D = np.load(paths[did])["Input"].astype(float); pp = sc11.prepare_sc(D, cm)
        sel, diag = qscan_w(D, cm, Q, pp)
        rec = dict(id=did, c_true=float(did.split("_")[1][1:]), lams=LAMS, qscan_w=sel, diag=diag,
                   seconds=time.monotonic() - t0)
    else:
        _, grp, well, mult, c, seed = sets[did]
        g = signal_L(ctx, W[well], c, seed) if grp == "L" else gen.signal(ctx, c, seed)
        D = g["x"] + mult * g["level"] * ctx["noise"][0]
        sr = sq.estimate(D, cm, Q, table=T, wins=wins, nboot=0, seed=seed)
        pp = sc11.prepare_sc(D, cm); fq, qi = sc11.qscan(D, cm, Q, pp)
        sel, diag = qscan_w(D, cm, Q, pp)
        rec = dict(id=did, group=grp, well=well, noise_mult=mult, c_true=c, seed=seed, level=g["level"], lines=ctx["lines"],
                   c_sr=float(sr["c"]), qscan=fq, qscan_boundary=bool(qi["at_boundary"]), lams=LAMS, qscan_w=sel,
                   diag=diag, seconds=time.monotonic() - t0)
    path = out / f"{did}.json"
    if not path.exists():
        path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
    (out / f"{did}.claim").unlink(missing_ok=True); print(did, f"{rec['seconds']:.0f}s", flush=True)
