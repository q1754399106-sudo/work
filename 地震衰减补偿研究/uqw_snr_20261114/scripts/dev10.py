"""Step 5 (development data, exploratory): UQW-S (uqws.py) with and without the target-band cap, all 24 development sets x
states of dev9.py, and FD H1-H6 (30, 14 dB). Methods UQW, HWin (= uqws.shape cap=False; equals dev9 HWin except above
fs_r), UQWS (cap=True), and VDBin (fair VDB: target w_in_r). usage: dev10.py semi ID [ID ...] | dev10.py fd -> results/dev10/"""
import sys, json, time
from env import np, semi, score, TW, snr_rows, CW, load, t2, lib12, SEL, tvfilt, CURVE, MS, H, DEV, ROOT, F10, noise_of, fd
import uqws as S

OUT = ROOT / "results" / "dev10"; OUT.mkdir(parents=True, exist_ok=True)


def methods(D, cm, Q, c, tw, anchor, w_t, fs_r, inv):
    X = CW.uqw(D, cm, Q, c, H.SIG, tw, anchor, H.GAMMA, inv_snr=inv)
    return dict(UQW=X, HWin=S.shape(X, w_t, inv, fs_r, cap=False), UQWS=S.shape(X, w_t, inv, fs_r),
                VDBin=lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, tw), w_t, SEL["B3"]["eps"], win=SEL["B3"]["win"]))


if sys.argv[1] == "semi":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"])
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); rec = []; sd = 20261114000 + DEV.index(did)
        for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
            for st, kind, x in (("res1", "res", 1.0), ("res4", "res", 4.0), ("gau4", "gau", 4.0)):
                D = sig + noise_of(arr, kind, x, sd); inv = snr_rows(D, TW, F10); w_t, fs_r = S.target(D, cm)
                for m in MS:
                    rec += [dict(state=var + st, method=k, m=m, S=score(Y, arr["ideal"])["S"], T2=t2(Y, arr["R"], fp))
                            for k, Y in methods(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, w_t, fs_r, inv).items()]
                rec.append(dict(state=var + st, method="_fs_r", m=0, fs_r=fs_r))
        (OUT / f"semi_{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
else:
    for h in fd.MODELS:
        D0, T, cm, tw, fp = fd.load_fd(h); rec = []; t0 = time.monotonic(); rk = fd.ricker(fp)
        for j, (lv, db) in enumerate((("30dB", 30.0), ("14dB", 14.0))):
            D = D0 + fd.noise(D0, db, 20261113900 + 10 * fd.MODELS.index(h) + j); inv = snr_rows(D, tw, F10)
            for m in MS:
                rec += [dict(state=lv, method=k, m=m, S=fd.score_fd(Y, T, tw)[0]) for k, Y in methods(D, cm, cm["Q"], m, tw, 0.0, rk, np.inf, inv).items()]
        (OUT / f"fd_{h}.json").write_text(json.dumps(rec), encoding="utf-8"); print(h, "%.0fs" % (time.monotonic() - t0), flush=True)
