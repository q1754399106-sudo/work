"""UQ-W confirmation (PLAN_CONFIRM_W_ZH.md): new semi-synthetic batch 5 (calibrated and 4x-noise versions) and new FD
models H7-H12 (14 and 30 dB). Operator frozen: uqw.uqw, gamma 0, sig 0.447, lognormal, 70 dB floor.

usage: confirmW.py run K N   worker K of N over 24 semi-synthetic + 36 FD jobs -> results/confirmW/<job>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uqw import np, semi, TW, uqw, snr_rows, DT, lib12, G_IQ  # noqa: E402
from uqiq import ROOT, BASE, uqiq  # noqa: E402
from stage0 import M_ALL, score, sc11, SEL  # noqa: E402
import fdcheck as fd  # noqa: E402
import scipy.io as sio  # noqa: E402

SIG, GAMMA, NOISE_X = 0.447, 0.0, 4.0
D5 = BASE / "semisynth_fault_q5_20261109"; FD2 = BASE / "fd2d_confirm_20261109"
IDS5 = [f"Z{k}_c{c}_r{r}" for k in "MNO" for c in (0.75, 1.05) for r in (0, 1)]
MODELS2 = [f"H{i}" for i in range(7, 13)]; LEVELS2 = {"14dB": 14.0, "30dB": 30.0}; R = (1, 2, 3)
OUT = ROOT / "results" / "confirmW"


def load5(did):
    meta = json.loads((D5 / "analysis" / "datasets.json").read_text(encoding="utf-8"))[did]
    return {k: v.astype(float) for k, v in np.load(D5 / "cache" / "data" / f"{did}.npz").items()}, meta


def load_fd2(name):
    m = sio.loadmat(FD2 / "models" / f"{name}.mat", squeeze_me=True, struct_as_record=False)
    t = sio.loadmat(FD2 / "truth" / f"{name}_truth.mat", squeeze_me=True, struct_as_record=False)
    c = m["config"]; cm = dict(Q=np.asarray(c.Q, float), vp=np.asarray(c.vp, float), h=np.asarray(c.h, float), fref=float(c.fref))
    return np.asarray(m["clean"], float), np.asarray(t["target"], float), cm, np.asarray(c.t, float), float(c.fpeak)


def methods(D, cm, Q, c, tw, anchor, w_out, inv):
    return dict(UQW=uqw(D, cm, Q, c, SIG, tw, anchor, GAMMA, inv_snr=inv),
                W0=uqw(D, cm, Q, c, 0.0, tw, anchor, 0.0, inv_snr=inv),
                VDB=lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, tw), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"]),
                IQ70=fd.iq_t(D, cm, Q, c, c, G_IQ, tw),
                UQIQ=fd.uq_t(D, cm, Q, c, SIG, tw, anchor),                    # descriptive
                C1=fd.iq_t(D, cm, Q, c, c, 30, tw))                              # descriptive


def run_semi(did, version, cm, Q):
    arr, meta = load5(did); D = arr["signal"] + (1.0 if version == "cal" else NOISE_X) * arr["noise"]; ideal = arr["ideal"]
    w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, np.fft.rfftfreq(1024, DT)); out = []
    for m in M_ALL:
        for k, X in methods(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, w_out, inv).items():
            out.append(dict(method=k, m=m, S=score(X, ideal)["S"]))
    return out


def run_fd(h, lv, r):
    D0, T, cm, tw, fp = load_fd2(h); i = MODELS2.index(h) + 1
    seed = 20261115000 + 1000 * i + 10 * list(LEVELS2).index(lv) + r
    D = D0 + fd.noise(D0, LEVELS2[lv], seed); inv = snr_rows(D, tw, np.fft.rfftfreq(1024, DT)); out = []
    for m in M_ALL:
        for k, X in methods(D, cm, cm["Q"], m, tw, 0.0, fd.ricker(fp), inv).items():
            out.append(dict(method=k, m=m, S=fd.score_fd(X, T, tw)[0]))
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
