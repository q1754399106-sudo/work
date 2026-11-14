"""Step 3b (development data, exploratory): UQW-Q (uqwq.py; QM mean shift only, QP posterior) in the states of dev7.py.
usage: dev8.py semi ID [ID ...] | dev8.py fd  -> results/dev8/   (UQW recomputed as the reference)
"""
import sys, json, time
from env import np, lib, semi, score, sc11, TW, DT, snr_rows, CW, load, t2, tvfilt, CURVE, MS, H, DEV, ROOT, F10, noise_of, fd
import uqwq as Qm

OUT = ROOT / "results" / "dev8"; OUT.mkdir(parents=True, exist_ok=True)


def three(D, cm, Q, c, tw, anchor, inv, A0, w0):
    return dict(UQW=CW.uqw(D, cm, Q, c, H.SIG, tw, anchor, H.GAMMA, inv_snr=inv), QM=Qm.uqwq(D, cm, Q, c, H.SIG, tw, anchor, inv, A0, w0, mode="mean"),
                QP=Qm.uqwq(D, cm, Q, c, H.SIG, tw, anchor, inv, A0, w0, mode="post"))


if sys.argv[1] == "semi":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"])
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); rec = []; sd = 20261114000 + DEV.index(did)
        for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
            for st, kind, x in (("res1", "res", 1.0), ("res4", "res", 4.0), ("gau4", "gau", 4.0)):
                D = sig + noise_of(arr, kind, x, sd); inv = snr_rows(D, TW, F10); A0, w0 = Qm.anchor_semi(D, TW)
                for m in MS:
                    rec += [dict(state=var + st, method=k, m=m, S=score(Y, arr["ideal"])["S"], T2=t2(Y, arr["R"], fp)) for k, Y in three(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, inv, A0, w0).items()]
        (OUT / f"semi_{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
else:
    for h in fd.MODELS:
        D0, T, cm, tw, fp = fd.load_fd(h); rec = []; t0 = time.monotonic(); A0 = np.abs(np.fft.rfft(fd.ricker(fp), Qm.NF)); w0 = np.ones_like(A0)
        for j, (lv, db) in enumerate((("30dB", 30.0), ("14dB", 14.0))):
            D = D0 + fd.noise(D0, db, 20261113900 + 10 * fd.MODELS.index(h) + j); inv = snr_rows(D, tw, F10)
            for m in MS:
                rec += [dict(state=lv, method=k, m=m, S=fd.score_fd(Y, T, tw)[0]) for k, Y in three(D, cm, cm["Q"], m, tw, 0.0, inv, A0, w0).items()]
        (OUT / f"fd_{h}.json").write_text(json.dumps(rec), encoding="utf-8"); print(h, "%.0fs" % (time.monotonic() - t0), flush=True)
