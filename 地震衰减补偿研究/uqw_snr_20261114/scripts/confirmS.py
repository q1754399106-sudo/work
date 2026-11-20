"""UQW-S confirmation (PLAN_CONFIRM_S_ZH.md): semi-synthetic batch 6 (calibrated, 4x residual, 4x spatially white noise; the
generator unchanged = variant O, plus the descriptive variant F) and FD H13-H18 (30 and 14 dB, 3 realisations).
Frozen code: uqws.py (UQW-S), uqw.py through confirmW (UQ-W: lognormal, sigma 0.447, gamma 0, 70 dB floor), lib12.tvw (VDB).
Methods: UQW, UQWS, VDBin (VDB, target = UQW-S's robust shallow wavelet), VDB (frozen target), HWin (UQW-S without the band
cap; descriptive), IQ70 (descriptive). Scores: semi S (lib15 score, 2.3-2.8 s) and T2 (descriptive); FD S (fdcheck.score_fd).
usage: confirmS.py smoke           code path and timing on one used development set and FD H1 (no scores printed)
       confirmS.py run K N         worker K of N over 12 semi + 36 FD jobs -> results/confirmS/<job>.json (resumable)
"""
import sys, json, time
from env import np, semi, score, TW, snr_rows, CW, t2, lib12, SEL, tvfilt, CURVE, H, ROOT, BASE, F10, fd
import uqws as S
import scipy.io as sio

MS = (0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0)
D6 = BASE / "semisynth_fault_q6_20261114"; FD3 = BASE / "fd2d_confirm3_20261114"
IDS6 = [f"Z{k}_c{c}_r{r}" for k in "PRT" for c in (0.75, 1.05) for r in (0, 1)]
MODELS3 = [f"H{i}" for i in range(13, 19)]; LEVELS3 = {"30dB": 30.0, "14dB": 14.0}; R3 = (1, 2, 3)
STATES = (("cal", "res", 1.0), ("res4", "res", 4.0), ("gau4", "gau", 4.0))
OUT = ROOT / "results" / "confirmS"


def gau_seed(i):
    return 20261180000 + i                                # i = position of the dataset in IDS6


def fd_seed(h, lv, r):
    return 20261190000 + 1000 * (MODELS3.index(h) + 1) + 10 * list(LEVELS3).index(lv) + r


def noise6(arr, kind, x, i):
    return x * (arr["noise"] if kind == "res" else semi.gaussian_noise(arr["noise"], np.random.default_rng(gau_seed(i))))


def methods(D, cm, Q, c, tw, anchor, w_t, fs_r, w_frozen, inv):
    X = CW.uqw(D, cm, Q, c, H.SIG, tw, anchor, H.GAMMA, inv_snr=inv); ph = fd.iq_t(D, cm, Q, None, c, 0, tw)
    v = lambda w: lib12.tvw(ph, w, SEL["B3"]["eps"], win=SEL["B3"]["win"])
    return dict(UQW=X, UQWS=S.shape(X, w_t, inv, fs_r), VDBin=v(w_t), VDB=v(w_frozen), HWin=S.shape(X, w_t, inv, fs_r, cap=False),
                IQ70=fd.iq_t(D, cm, Q, c, c, CW.G_IQ, tw))


def run_semi(arr, meta, i, cm, Q, quiet=False):
    fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"]); rec, bands = [], []
    for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
        for st, kind, x in STATES:
            D = sig + noise6(arr, kind, x, i); inv = snr_rows(D, TW, F10); w_t, fs_r = S.target(D, cm); wv = CW.sc11.prepare_sc(D, cm)["wv"]
            bands.append(dict(var=var, state=st, fs_r=fs_r, fs=wv["fmax_shallow"], fd=wv["fmax_target"]))
            for m in MS:
                for k, Y in methods(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, w_t, fs_r, wv["w_out"], inv).items():
                    rec.append(dict(var=var, state=st, method=k, m=m, S=score(Y, arr["ideal"])["S"], T2=t2(Y, arr["R"], fp), finite=bool(np.isfinite(Y).all())))
    return rec, bands


def load_fd3(name):
    m = sio.loadmat(FD3 / "models" / f"{name}.mat", squeeze_me=True, struct_as_record=False)
    t = sio.loadmat(FD3 / "truth" / f"{name}_truth.mat", squeeze_me=True, struct_as_record=False)
    c = m["config"]; cm = dict(Q=np.asarray(c.Q, float), vp=np.asarray(c.vp, float), h=np.asarray(c.h, float), fref=float(c.fref))
    return np.asarray(m["clean"], float), np.asarray(t["target"], float), cm, np.asarray(c.t, float), float(c.fpeak)


def run_fd_unit(D0, T, cm, tw, fp, db, seed):
    D = D0 + fd.noise(D0, db, seed); inv = snr_rows(D, tw, F10); rk = fd.ricker(fp); rec = []
    for m in MS:
        for k, Y in methods(D, cm, cm["Q"], m, tw, 0.0, rk, np.inf, rk, inv).items():
            rec.append(dict(method=k, m=m, S=fd.score_fd(Y, T, tw)[0], finite=bool(np.isfinite(Y).all())))
    return rec


if __name__ == "__main__" and sys.argv[1] == "smoke":
    from env import load
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); t0 = time.monotonic()
    arr, meta = load("ZA_c0.75_r0"); rec, bands = run_semi(arr, meta, 0, cm, Q)
    print("semi smoke: %d records, all finite %s, states %s, %.0fs" % (len(rec), all(r["finite"] for r in rec), sorted({(r['var'], r['state']) for r in rec}), time.monotonic() - t0))
    print("D == signal + noise: max |d| %.1e" % np.abs(arr["D"] - (arr["signal"] + arr["noise"])).max())
    t0 = time.monotonic(); D0, T, cmf, tw, fp = fd.load_fd("H1"); rec = run_fd_unit(D0, T, cmf, tw, fp, 30.0, 1)
    print("FD smoke: %d records, all finite %s, %.0fs" % (len(rec), all(r["finite"] for r in rec), time.monotonic() - t0))
    print("seeds: gau", gau_seed(0), "..", gau_seed(11), " fd", fd_seed("H13", "30dB", 1), "..", fd_seed("H18", "14dB", 3))

if __name__ == "__main__" and sys.argv[1] == "run":
    sys.path.insert(0, str(D6 / "scripts")); import fq6  # noqa: E402
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [("semi", d) for d in IDS6] + [("fd", (h, lv, r)) for h in MODELS3 for lv in LEVELS3 for r in R3]
    for j, (kind, key) in enumerate(jobs):
        name = "S6_" + key if kind == "semi" else "FD3_%s_%s_r%d" % key; path = OUT / f"{name}.json"
        if j % n != k or path.exists():
            continue
        t0 = time.monotonic()
        if kind == "semi":
            arr, meta = fq6.load(key); arr = {a: v.astype(float) for a, v in arr.items()}
            recs, bands = run_semi(arr, meta, IDS6.index(key), cm, Q); out = dict(id=name, records=recs, bands=bands)
        else:
            h, lv, r = key; D0, T, cmf, tw, fp = load_fd3(h); s = fd_seed(h, lv, r); out = dict(id=name, seed=s, records=run_fd_unit(D0, T, cmf, tw, fp, LEVELS3[lv], s))
        out["seconds"] = time.monotonic() - t0; path.write_text(json.dumps(out), encoding="utf-8"); print(name, len(out["records"]), "%.0fs" % out["seconds"], flush=True)
