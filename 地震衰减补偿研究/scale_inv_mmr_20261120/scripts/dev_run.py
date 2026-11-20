"""Development run (PLAN_ZH.md): the stage-W development data of qunknown_comp_20261109 (semi-synthetic DEV, 24 datasets;
FD H1-H6 at 14 dB, 3 realisations, stage 2b seeds). For each Q multiplier m: the new rules SIU, SIR (si.py) and the
baselines recomputed so that T2 can be scored - UQW (uqw.uqw, gamma 0), W0 (sigma 0), VDB, UQWl / CBL0.5 (cbl.design_cbl)
and RGJ (mmr operator rebuilt from the archived prior of results/mmr; mmr.design costs about 140 s per m). Stored S of UQW, CBL0.5 and RGJ in qunknown_comp results are re-checked by analyze_dev.py.
Scores: S, M1, M1_lag (both data types); T2 (semi-synthetic only, probe.t2).

usage: dev_run.py run K N   worker K of N over 24 + 18 jobs -> results/dev/<job>.json
"""
import sys, json, time
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import si  # noqa: E402
from si import np, semi, DT, NFFT, snr_rows, MAIN  # noqa: E402
from uqiq import TW, load  # noqa: E402
from stage0 import DEV, M_ALL, score, SEL, lib  # noqa: E402
from uqw import uqw  # noqa: E402
import mmr, cbl  # noqa: E402
import confirmW as CW  # noqa: E402
import fdcheck as fd  # noqa: E402
import lib12  # noqa: E402

OUT = ROOT / "results" / "dev"; F = np.fft.rfftfreq(NFFT, DT); SIG = 0.447
WS = semi.win(TW, 2.3, 2.8)


def t2(X, R, fp):
    return float(np.mean([lib.tie(X[:, j], semi.ricker_ref(R[:, j], fp), TW, WS)[0] for j in range(2, X.shape[1], 5)]))


def rules(D, cm, Q, c, tw, anchor, inv, w, w_out, phase_iq, rgj_pi=None):
    X, info = {}, {}
    X["UQW"] = uqw(D, cm, Q, c, SIG, tw, anchor, 0.0, inv_snr=inv)
    X["W0"] = uqw(D, cm, Q, c, 0.0, tw, anchor, 0.0, inv_snr=inv)
    X["VDB"] = lib12.tvw(phase_iq(D, cm, Q, c, tw), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"])
    outs, ci = cbl.design_cbl(D, cm, Q, c, tw, anchor, inv, taus=(1.0, 0.5)); X["UQWl"], X["CBL0.5"] = outs[1.0], outs[0.5]
    info["CBL"] = dict(R_mm=ci["R_mm"], R_T=ci["R_T"], lam=np.round(ci["runs"]["0.5"]["lam"], 5).tolist())
    if rgj_pi is not None:                                    # RGJ rebuilt from the archived least-favourable prior (S re-checked)
        X["RGJ"] = rgj(D, cm, Q, c, tw, anchor, inv, rgj_pi)
    for r in ("SIU", "SIR"):
        t0 = time.monotonic(); X[r], inf = si.design(D, cm, Q, c, tw, anchor, r, inv, w); inf["sec"] = time.monotonic() - t0
        info[r] = {k: (np.round(v, 4).tolist() if isinstance(v, list) else v) for k, v in inf.items()}
    return X, info


QC = ROOT.parent / "qunknown_comp_20261109" / "results" / "mmr"           # tracked results (this checkout)


def archived_pi(name):
    r = json.loads((QC / f"{name}.json").read_text(encoding="utf-8"))
    return {x["m"]: np.asarray(x["info"]["pi"], float) for x in r["records"] if x["method"] == "RGJ"}


def rgj(D, cm, Q, c, tw, anchor, inv, pi):
    """mmr.design's RGJ operator for a given prior (renormalised after the 4-decimal rounding of the archive)."""
    f, lf, x_hat, dx_hat = mmr.setup(cm, Q, c, tw, anchor); s = semi.gain_sigma2(mmr.G_IQ)
    xs, _ = mmr.profiles(x_hat, dx_hat, "scalar"); Gc, den, _, sc = mmr.cells(xs, f, lf, s, inv)
    return mmr.apply_rows(mmr.h_of(pi / pi.sum(), Gc, den) * sc, D)


def run_semi(did, cm, Q):
    arr, meta = load(did); D, ideal, R = arr["D"], arr["ideal"], arr["R"]; inv = snr_rows(D, TW, F); fp = float(cm["fpeak"])
    w_out = CW.sc11.prepare_sc(D, cm)["wv"]["w_out"]; w = si.wspec(w_out); out = []
    piq = lambda D, cm, Q, c, tw: lib12.iq(D, cm, Q, None, c, 0); P = archived_pi(did)
    for m in M_ALL:
        X, info = rules(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, inv, w, w_out, piq, P[m])
        for k, v in X.items():
            s = score(v, ideal); out.append(dict(method=k, m=m, S=s["S"], M1=s["M1"], M1_lag=s["M1_lag"], T2=t2(v, R, fp)))
        out.append(dict(method="info", m=m, info=info))
    return out


def run_fd(h, r):
    D0, T, cm, tw, fp = fd.load_fd(h); i = fd.MODELS.index(h) + 1
    seed = 20261113000 + 1000 * i + 10 * list(fd.LEVELS).index("14dB") + r          # stage 2b / W seeds
    D = D0 + fd.noise(D0, fd.LEVELS["14dB"], seed); inv = snr_rows(D, tw, F); w_out = fd.ricker(fp); w = si.wspec(w_out); out = []
    piq = lambda D, cm, Q, c, tw: fd.iq_t(D, cm, Q, None, c, 0, tw); P = archived_pi("FD_%s_14dB_r%d" % (h, r))
    for m in M_ALL:
        X, info = rules(D, cm, cm["Q"], m, tw, 0.0, inv, w, w_out, piq, P[m])
        for k, v in X.items():
            S, z, lag = fd.score_fd(v, T, tw); out.append(dict(method=k, m=m, S=S, M1=z, M1_lag=lag))
        out.append(dict(method="info", m=m, info=info))
    return out, seed


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
            R_, seed = run_fd(*key); rec = dict(id=name, seed=seed, records=R_)
        rec["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(name, f"{rec['seconds']:.0f}s", flush=True)
