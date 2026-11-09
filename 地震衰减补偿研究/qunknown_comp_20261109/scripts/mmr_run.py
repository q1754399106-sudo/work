"""Minimax-regret prototype (PLAN_MMR_ZH.md) on the stage-W development data: semi-synthetic DEV (24 datasets) and
FD H1-H6 at 14 dB (3 realisations, stage 2b seeds). Rules UQW, MMJ, RGC, RGJ, RGJ2 from mmr.py; VDB and W0 are
taken from stages 0, 2b and W (identical data and noise).

usage: mmr_run.py run K N   worker K of N over 24 + 18 jobs -> results/mmr/<job>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mmr import np, semi, DT, NFFT, snr_rows, design  # noqa: E402
from uqiq import ROOT, TW, load  # noqa: E402
from stage0 import DEV, M_ALL, score  # noqa: E402
import fdcheck as fd  # noqa: E402

RULES = ("UQW", "MMJ", "RGJ", "RGC", "RGJ2"); OUT = ROOT / "results" / "mmr"


def slim(info):
    return {k: (np.round(v, 4).tolist() if isinstance(v, list) else v) for k, v in info.items()}


def run_semi(did, cm, Q):
    arr, meta = load(did); D, ideal = arr["D"], arr["ideal"]; inv = snr_rows(D, TW, np.fft.rfftfreq(NFFT, DT)); R = []
    for m in M_ALL:
        for rule in RULES:
            X, info = design(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, rule, inv)
            R.append(dict(method=rule, m=m, S=score(X, ideal)["S"], info=slim(info)))
    return R


def run_fd(h, r):
    D0, T, cm, tw, fp = fd.load_fd(h); i = fd.MODELS.index(h) + 1
    seed = 20261113000 + 1000 * i + 10 * list(fd.LEVELS).index("14dB") + r          # stage 2b / W seeds
    D = D0 + fd.noise(D0, fd.LEVELS["14dB"], seed); inv = snr_rows(D, tw, np.fft.rfftfreq(NFFT, DT)); R = []
    for m in M_ALL:
        for rule in RULES:
            X, info = design(D, cm, cm["Q"], m, tw, 0.0, rule, inv)
            R.append(dict(method=rule, m=m, S=fd.score_fd(X, T, tw)[0], info=slim(info)))
    return R, seed


if __name__ == "__main__" and sys.argv[1] == "run":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [("semi", d) for d in DEV] + [("fd", (h, r)) for h in fd.MODELS for r in fd.R]
    for j, (kind, key) in enumerate(jobs):
        name = key if kind == "semi" else "FD_%s_14dB_r%d" % key; path = OUT / f"{name}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic()
        if kind == "semi":
            rec = dict(id=name, records=run_semi(key, cm, Q))
        else:
            R, seed = run_fd(*key); rec = dict(id=name, seed=seed, records=R)
        rec["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(name, len(rec["records"]), f"{rec['seconds']:.0f}s", flush=True)
