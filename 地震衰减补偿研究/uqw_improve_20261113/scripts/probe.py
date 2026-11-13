"""Step 3 probe (development data only): do the used semi-synthetic development sets reproduce the field behaviour?

For each dataset and Q multiplier: S (score against the ideal section, as stage 0), a field-style tie T2 (max correlation
within +-20 ms of the output trace with the true reflectivity * Ricker(fpeak), window 2.3-2.8 s, every 5th trace; lib.tie and
semi.ricker_ref as on the field wells) and the share of output energy above 80 Hz, for UQ-W, VDB and IQ70 (confirmW.methods).

usage: probe.py ID [ID ...]
"""
import sys, json, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "qunknown_comp_20261109" / "scripts"))
import confirmW as CW  # noqa: E402
from confirmW import np, semi, score, sc11, TW, DT, snr_rows  # noqa: E402
from stage0 import load, lib  # noqa: E402

MS = (0.5, 0.7, 1.0, 1.4, 2.0); KEYS = ("UQW", "VDB", "IQ70"); F = np.fft.rfftfreq(1024, DT); WS = semi.win(TW, 2.3, 2.8)


def t2(X, R, fp):
    return float(np.mean([lib.tie(X[:, j], semi.ricker_ref(R[:, j], fp), TW, WS)[0] for j in range(2, X.shape[1], 5)]))


def hf(X):
    k = (TW >= 2.0) & (TW <= 2.8); P = np.mean(np.abs(np.fft.rfft(X[k] * np.hanning(k.sum())[:, None], 1024, axis=0)) ** 2, axis=1); return float(P[F > 80].sum() / P.sum())


if __name__ == "__main__":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); out = ROOT / "results" / "probe"; out.mkdir(parents=True, exist_ok=True)
    for did in sys.argv[1:]:
        t0 = time.monotonic(); arr, meta = load(did); D, ideal, R = arr["D"], arr["ideal"], arr["R"]; w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, F)
        rec = [dict(method="Input", m=1.0, S=score(D, ideal)["S"], T2=t2(D, R, fp), hf=hf(D)), dict(method="Ideal", m=1.0, S=1.0, T2=t2(ideal, R, fp), hf=hf(ideal))]
        for m in MS:
            r = CW.methods(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, w_out, inv)
            rec += [dict(method=k, m=m, S=score(r[k], ideal)["S"], T2=t2(r[k], R, fp), hf=hf(r[k])) for k in KEYS]
        (out / f"{did}.json").write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, "fpeak %.1f  %.0fs" % (fp, time.monotonic() - t0))
        for q in rec[:2]:
            print("  %-5s            S %.3f  T2 %.3f  >80Hz %.3f" % (q["method"], q["S"], q["T2"], q["hf"]))
        for k in KEYS:
            for name in ("S", "T2", "hf"):
                print("  %-5s %-3s by m %s: " % (k, name, MS) + " ".join("%.3f" % q[name] for q in rec if q["method"] == k), flush=True)
