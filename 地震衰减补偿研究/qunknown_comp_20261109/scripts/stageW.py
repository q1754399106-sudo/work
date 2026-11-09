"""Stage W (PLAN_STAGEW_ZH.md): UQ-W (data-driven noise term) on the semi-synthetic development set and the FD data.

usage: stageW.py run K N   worker K of N over 24 semi-synthetic + 36 FD jobs -> results/stageW/<job>.json
Other methods (IQ70, VDB, UQ-IQ, C1, C2) are taken from stages 0, 1 and 2b (identical data and noise).
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqw import np, semi, TW, uqw, snr_rows, DT  # noqa: E402
from uqiq import ROOT, load  # noqa: E402
from stage0 import DEV, M_ALL, score  # noqa: E402
import fdcheck as fd  # noqa: E402

SIG = 0.447; SIG_DESC = (0.245, 0.894); OUT = ROOT / "results" / "stageW"
VARIANTS = [("W0", 0.0, 0.0)] + [("UQW_g0", SIG, 0.0), ("UQW_g1", SIG, 1.0)] + [("UQW_g0", q, 0.0) for q in SIG_DESC]


def run_semi(did, cm, Q):
    arr, meta = load(did); D, ideal = arr["D"], arr["ideal"]; ce = meta["c_eff_cum"]
    inv = snr_rows(D, TW, np.fft.rfftfreq(1024, DT)); R = []
    for m in M_ALL:
        for name, sig, g in VARIANTS:
            R.append(dict(method=name, m=m, sig=sig, gamma=g, S=score(uqw(D, cm, Q, m * ce, sig, TW, 1.4, g, inv_snr=inv), ideal)["S"]))
    return R


def run_fd(h, lv, r):
    D0, T, cm, tw, fp = fd.load_fd(h); i = fd.MODELS.index(h) + 1
    seed = 20261113000 + 1000 * i + 10 * list(fd.LEVELS).index(lv) + r                      # as stage 2b
    D = D0 + fd.noise(D0, fd.LEVELS[lv], seed); inv = snr_rows(D, tw, np.fft.rfftfreq(1024, DT)); R = []
    for m in M_ALL:
        for name, sig, g in VARIANTS:
            R.append(dict(method=name, m=m, sig=sig, gamma=g, S=fd.score_fd(uqw(D, cm, cm["Q"], m, sig, tw, 0.0, g, inv_snr=inv), T, tw)[0]))
    return R, seed


if __name__ == "__main__" and sys.argv[1] == "run":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [("semi", d) for d in DEV] + [("fd", (h, lv, r)) for h in fd.MODELS for lv in fd.LEVELS for r in fd.R]
    for j, (kind, key) in enumerate(jobs):
        name = key if kind == "semi" else "FD_%s_%s_r%d" % key; path = OUT / f"{name}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic()
        if kind == "semi":
            rec = dict(id=key, records=run_semi(key, cm, Q))
        else:
            R, seed = run_fd(*key); rec = dict(id=name, seed=seed, records=R)
        rec["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(name, len(rec["records"]), f"{rec['seconds']:.0f}s", flush=True)
