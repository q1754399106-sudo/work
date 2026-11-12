"""CBL confirmation (PLAN_CONFIRM_CBL_ZH.md): CBL frozen (cbl.design_cbl, 8b5836a) on the UQ-W confirmation data, which CB/CBL
never touched: semi-synthetic batch 5 (calibrated and 4x noise) and FD H7-H12 (14 and 30 dB, confirmW.py seeds).
Rules UQWl, CBL0.5 (confirmed), CBL0.25 and UQW_re (uqw.uqw recomputed; must equal results/confirmW UQW).
UQW, VDB, IQ70 and W0 are read from results/confirmW by analyze_confirm_cbl.py.

usage: confirm_cbl.py run K N   worker K of N over 24 semi-synthetic + 36 FD jobs -> results/confirm_cbl/<job>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cbl import np, semi, DT, NFFT, snr_rows, design_cbl  # noqa: E402
from uqw import uqw, TW  # noqa: E402
from uqiq import ROOT  # noqa: E402
from stage0 import M_ALL, score  # noqa: E402
from confirmW import load5, load_fd2, IDS5, MODELS2, LEVELS2, R, NOISE_X, SIG  # noqa: E402
import fdcheck as fd  # noqa: E402

OUT = ROOT / "results" / "confirm_cbl"; NAME = {1.0: "UQWl", 0.5: "CBL0.5", 0.25: "CBL0.25"}


def slim(info):
    return {k: (slim(v) if isinstance(v, dict) else np.round(v, 5).tolist() if isinstance(v, list) else v) for k, v in info.items()}


def rules(D, cm, Q, c, tw, anchor, inv):
    outs, info = design_cbl(D, cm, Q, c, tw, anchor, inv); X = {NAME[t]: v for t, v in outs.items()}
    X["UQW_re"] = uqw(D, cm, Q, c, SIG, tw, anchor, 0.0, inv_snr=inv); return X, info


def run_semi(did, version, cm, Q):
    arr, meta = load5(did); D = arr["signal"] + (1.0 if version == "cal" else NOISE_X) * arr["noise"]; ideal = arr["ideal"]
    inv = snr_rows(D, TW, np.fft.rfftfreq(NFFT, DT)); out = []
    for m in M_ALL:
        X, info = rules(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, inv)
        out += [dict(method=k, m=m, S=score(v, ideal)["S"]) for k, v in X.items()] + [dict(method="info", m=m, info=slim(info))]
    return out


def run_fd(h, lv, r):
    D0, T, cm, tw, fp = load_fd2(h); i = MODELS2.index(h) + 1
    seed = 20261115000 + 1000 * i + 10 * list(LEVELS2).index(lv) + r                  # confirmW.py seeds
    D = D0 + fd.noise(D0, LEVELS2[lv], seed); inv = snr_rows(D, tw, np.fft.rfftfreq(NFFT, DT)); out = []
    for m in M_ALL:
        X, info = rules(D, cm, cm["Q"], m, tw, 0.0, inv)
        out += [dict(method=k, m=m, S=fd.score_fd(v, T, tw)[0]) for k, v in X.items()] + [dict(method="info", m=m, info=slim(info))]
    return out, seed


if __name__ == "__main__" and sys.argv[1] == "run":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [("semi", (d, v)) for d in IDS5 for v in ("cal", "nd")] + [("fd", (h, lv, r)) for h in MODELS2 for lv in LEVELS2 for r in R]
    for j, (kind, key) in enumerate(jobs):
        name = "S5_%s_%s" % key if kind == "semi" else "FD2_%s_%s_r%d" % key; path = OUT / f"{name}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic()
        if kind == "semi":
            rec = dict(id=name, records=run_semi(*key, cm, Q))
        else:
            recs, seed = run_fd(*key); rec = dict(id=name, seed=seed, records=recs)
        rec["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(name, len(rec["records"]), f"{rec['seconds']:.0f}s", flush=True)
