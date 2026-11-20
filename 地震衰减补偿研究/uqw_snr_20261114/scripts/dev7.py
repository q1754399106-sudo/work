"""Step 3 (development data, exploratory): UQW-P (uqwp.py) against UQ-W, UQW-H and VDB in every development state.

semi: variants O and F (probe2.tvfilt), states res1 (calibrated), res4, gau4 (4x noise, real residual / spatially white).
  UQW, VDB, H (uqwh.half), P50 (const 0.5), PR (reliability weights), PR1 (lam = w w0), and with the oracle clean target
  (wavelet of the 1x data; only informative at 4x) H_tc, VDB_tc.
fd: H1-H6 at 30 and 14 dB (seeds of dev4.py); anchor t = 0 with the known Ricker spectrum, w0 = 1.
usage: dev7.py semi ID [ID ...] | dev7.py fd  -> results/dev7/
"""
import sys, json, time
from env import np, lib, semi, score, sc11, TW, DT, snr_rows, CW, load, t2, lib12, SEL, tvfilt, CURVE, MS, H, DEV, ROOT, F10, noise_of, fd
import uqwp as U

OUT = ROOT / "results" / "dev7"; OUT.mkdir(parents=True, exist_ok=True)
VC = np.asarray(json.loads((ROOT / "results" / "vcolour.json").read_text(encoding="utf-8"))["vc"])


def cands(D, X, cm, Q, c, tw, anchor, inv, A0, w0):
    dx = U.dxhat(cm, Q, c, tw, anchor); k = dict(sig=H.SIG, vc=VC)
    return dict(P50=U.update(X, D, tw, dx, inv, A0, w0, mode="const", **k), PR=U.update(X, D, tw, dx, inv, A0, w0, mode="rel", **k),
                PR1=U.update(X, D, tw, dx, inv, A0, w0, mode="w", **k))


if sys.argv[1] == "semi":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"])
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); rec = []; sd = 20261114000 + DEV.index(did)
        for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
            w_cl = sc11.prepare_sc(sig + arr["noise"], cm)["wv"]["w_out"]
            for st, kind, x in (("res1", "res", 1.0), ("res4", "res", 4.0), ("gau4", "gau", 4.0)):
                D = sig + noise_of(arr, kind, x, sd); w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, F10); A0, w0 = U.anchor_semi(D, TW)
                for m in MS:
                    c = m * meta["c_eff_cum"]; X = CW.uqw(D, cm, Q, c, H.SIG, TW, 1.4, H.GAMMA, inv_snr=inv); ph = fd.iq_t(D, cm, Q, None, c, 0, TW)
                    o = dict(UQW=X, VDB=lib12.tvw(ph, w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"]), H=H.half(X, w_out), **cands(D, X, cm, Q, c, TW, 1.4, inv, A0, w0))
                    if x > 1:
                        o.update(H_tc=H.half(X, w_cl), VDB_tc=lib12.tvw(ph, w_cl, SEL["B3"]["eps"], win=SEL["B3"]["win"]))
                    rec += [dict(state=var + st, method=k, m=m, S=score(Y, arr["ideal"])["S"], T2=t2(Y, arr["R"], fp)) for k, Y in o.items()]
        (OUT / f"semi_{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
else:
    f5 = np.fft.rfftfreq(U.NF, DT)
    for h in fd.MODELS:
        D0, T, cm, tw, fp = fd.load_fd(h); rec = []; t0 = time.monotonic(); A0 = np.abs(np.fft.rfft(fd.ricker(fp), U.NF)); w0 = np.ones_like(A0)
        for j, (lv, db) in enumerate((("30dB", 30.0), ("14dB", 14.0))):
            D = D0 + fd.noise(D0, db, 20261113900 + 10 * fd.MODELS.index(h) + j); inv = snr_rows(D, tw, F10)
            for m in MS:
                X = CW.uqw(D, cm, cm["Q"], m, H.SIG, tw, 0.0, H.GAMMA, inv_snr=inv)
                o = dict(UQW=X, VDB=lib12.tvw(fd.iq_t(D, cm, cm["Q"], None, m, 0, tw), fd.ricker(fp), SEL["B3"]["eps"], win=SEL["B3"]["win"]), H=H.half(X, fd.ricker(fp)),
                         **cands(D, X, cm, cm["Q"], m, tw, 0.0, inv, A0, w0))
                rec += [dict(state=lv, method=k, m=m, S=fd.score_fd(Y, T, tw)[0]) for k, Y in o.items()]
        (OUT / f"fd_{h}.json").write_text(json.dumps(rec), encoding="utf-8"); print(h, "%.0fs" % (time.monotonic() - t0), flush=True)
