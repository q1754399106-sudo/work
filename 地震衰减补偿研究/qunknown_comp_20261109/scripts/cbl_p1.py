"""Pre-freeze development check P1 of PLAN_CONFIRM_CBL_ZH.md: CBL on the two noise states the CB/CBL development never saw,
on development data only:
  semi-synthetic DEV (24) with 4x noise (D = signal + 4 noise, ideal of the calibrated version), and
  FD H1-H6 at 30 dB, 3 realisations, seed 20261113000 + 1000 i + 10 * 2 + r (stage 2b formula, level index 2 after 14dB, 08dB).
Rules: UQWl, CBL0.5, CBL0.25 (cbl.design_cbl, frozen) and the controls UQW (uqw.uqw, gamma 0) and VDB (as confirmW.py).
Used only to write the predictions of the confirmation plan; thresholds are fixed before this run.

usage: cbl_p1.py run K N   worker K of N over 24 + 18 jobs -> results/cbl_p1/<job>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cbl import np, semi, DT, NFFT, snr_rows, design_cbl  # noqa: E402
from uqw import uqw, lib12  # noqa: E402
from uqiq import ROOT, TW, load  # noqa: E402
from stage0 import DEV, M_ALL, score, sc11, SEL  # noqa: E402
import fdcheck as fd  # noqa: E402

OUT = ROOT / "results" / "cbl_p1"; NAME = {1.0: "UQWl", 0.5: "CBL0.5", 0.25: "CBL0.25"}; SIG, NOISE_X, DB30 = 0.447, 4.0, 30.0


def slim(info):
    return {k: (slim(v) if isinstance(v, dict) else np.round(v, 5).tolist() if isinstance(v, list) else v) for k, v in info.items()}


def rules(D, cm, Q, c, tw, anchor, inv, w_out):
    outs, info = design_cbl(D, cm, Q, c, tw, anchor, inv)
    X = {NAME[t]: v for t, v in outs.items()}
    X["UQW"] = uqw(D, cm, Q, c, SIG, tw, anchor, 0.0, inv_snr=inv)
    X["VDB"] = lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, tw), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"])
    return X, info


def run_semi(did, cm, Q):
    arr, meta = load(did); D = arr["signal"] + NOISE_X * arr["noise"]; ideal = arr["ideal"]
    w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, np.fft.rfftfreq(NFFT, DT)); R = []
    for m in M_ALL:
        X, info = rules(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, inv, w_out)
        R += [dict(method=k, m=m, S=score(v, ideal)["S"]) for k, v in X.items()] + [dict(method="info", m=m, info=slim(info))]
    return R


def run_fd(h, r):
    D0, T, cm, tw, fp = fd.load_fd(h); i = fd.MODELS.index(h) + 1; seed = 20261113000 + 1000 * i + 10 * 2 + r
    D = D0 + fd.noise(D0, DB30, seed); inv = snr_rows(D, tw, np.fft.rfftfreq(NFFT, DT)); R = []
    for m in M_ALL:
        X, info = rules(D, cm, cm["Q"], m, tw, 0.0, inv, fd.ricker(fp))
        R += [dict(method=k, m=m, S=fd.score_fd(v, T, tw)[0]) for k, v in X.items()] + [dict(method="info", m=m, info=slim(info))]
    return R, seed


if __name__ == "__main__" and sys.argv[1] == "run":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [("semi", d) for d in DEV] + [("fd", (h, r)) for h in fd.MODELS for r in fd.R]
    for j, (kind, key) in enumerate(jobs):
        name = f"{key}_x4" if kind == "semi" else "FD_%s_30dB_r%d" % key; path = OUT / f"{name}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic()
        if kind == "semi":
            rec = dict(id=name, records=run_semi(key, cm, Q))
        else:
            R, seed = run_fd(*key); rec = dict(id=name, seed=seed, records=R)
        rec["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(name, len(rec["records"]), f"{rec['seconds']:.0f}s", flush=True)
