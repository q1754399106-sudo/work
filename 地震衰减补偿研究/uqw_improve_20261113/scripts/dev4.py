"""Step 3, fourth probe (development data only): UQW-H (uqwh.py) in the states that the earlier probes did not cover.
  nd ID...  noise-dominated versions (signal + 4 noise) of the development sets, variants O and F: S and T2
  fd        finite-difference development models H1-H6 (fd2d_20261025), 30 and 14 dB, realisation 1: S (fdcheck.score_fd)
Methods: UQW, VDB, UQWH. usage: dev4.py nd ID [ID ...] | dev4.py fd
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe2 import np, semi, score, sc11, TW, DT, snr_rows, CW, load, ROOT, t2, lib12, SEL, tvfilt, CURVE, MS  # noqa: E402
import uqwh as H  # noqa: E402
fd = CW.fd; OUT = ROOT / "results" / "dev4"; OUT.mkdir(parents=True, exist_ok=True); F = np.fft.rfftfreq(1024, DT)


def three(D, cm, Q, c, tw, anchor, w_out):
    inv = snr_rows(D, tw, F); U = CW.uqw(D, cm, Q, c, H.SIG, tw, anchor, H.GAMMA, inv_snr=inv)
    return dict(UQW=U, VDB=lib12.tvw(fd.iq_t(D, cm, Q, None, c, 0, tw), w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"]), UQWH=H.half(U, w_out))


if sys.argv[1] == "nd":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"])
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); rec = []
        for var, sig in (("O", arr["signal"]), ("F", tvfilt(arr["signal"], r))):
            D = sig + 4.0 * arr["noise"]; w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]
            for m in MS:
                rec += [dict(var=var, method=k, m=m, S=score(X, arr["ideal"])["S"], T2=t2(X, arr["R"], fp)) for k, X in three(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, w_out).items()]
        (OUT / f"nd_{did}.json").write_text(json.dumps(rec), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
else:
    for h in fd.MODELS:
        D0, T, cm, tw, fp = fd.load_fd(h); rec = []; t0 = time.monotonic()
        for j, (lv, db) in enumerate((("30dB", 30.0), ("14dB", 14.0))):
            D = D0 + fd.noise(D0, db, 20261113900 + 10 * fd.MODELS.index(h) + j)
            for m in MS:
                rec += [dict(level=lv, method=k, m=m, S=fd.score_fd(X, T, tw)[0]) for k, X in three(D, cm, cm["Q"], m, tw, 0.0, fd.ricker(fp)).items()]
        (OUT / f"fd_{h}.json").write_text(json.dumps(rec), encoding="utf-8"); print(h, "fpeak %.1f n %d %.0fs" % (fp, len(tw), time.monotonic() - t0), flush=True)
