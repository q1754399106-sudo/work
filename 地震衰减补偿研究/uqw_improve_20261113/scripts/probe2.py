"""Step 3, second probe (development data only): a development variant whose depth-dependent spectral change is not of
constant-Q form, and a first look at UQ-W followed by the VDB shaping.

Variant "F": signal' = time-varying zero-phase filter of the signal with log gain r(f) (t - 1.4) / 1.15, where r(f) is the
departure of the field log spectral ratio (2.3-2.8 s over 1.2-1.6 s, 15 Volve lines, no wells) from its straight-line fit over
10-55 Hz (specform.py; 11-bin smoothing, held constant below 8 and above 72 Hz). D' = signal' + noise; ideal and R unchanged.
Variant "O": the original development data. Methods: UQ-W, VDB, IQ70 (confirmW.methods) and UQWS = lib12.tvw(UQ-W) with the
VDB parameters. Scores as probe.py.

usage: probe2.py curve            measure r(f) on the field lines -> results/field_residual.json
       probe2.py run ID [ID ...]  -> results/probe2/<ID>.json
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import probe as P  # noqa: E402
from probe import np, semi, score, sc11, TW, DT, snr_rows, CW, load, lib, ROOT, BASE, t2, hf  # noqa: E402
from confirmW import lib12, SEL  # noqa: E402

NF = 512; F5 = np.fft.rfftfreq(NF, DT); CURVE = ROOT / "results" / "field_residual.json"; MS = P.MS; KEYS = ("UQW", "VDB", "IQ70", "UQWS")


def tvfilt(X, r, win=0.4, hop=0.1):
    n = X.shape[0]; Lw = int(round(win / DT)); step = int(round(hop / DT)); h = np.hanning(Lw); out = np.zeros_like(X); ws = np.zeros(n)
    for s0 in range(-Lw // 2, n - Lw // 2, step):
        i0, i1 = max(s0, 0), min(s0 + Lw, n); hw = h[i0 - s0:i1 - s0]; tc = TW[0] + DT * (s0 + Lw / 2); g = np.exp(r * (tc - 1.4) / 1.15)
        out[i0:i1] += np.fft.irfft(np.fft.rfft(X[i0:i1] * hw[:, None], NF, axis=0) * g[:, None], NF, axis=0)[:i1 - i0]; ws[i0:i1] += hw
    return out / np.maximum(ws, 1e-3 * ws.max())[:, None]


if __name__ == "__main__" and sys.argv[1] == "curve":
    import subprocess
    E13 = BASE / "field_baselines_volve_20261029"; sys.path.insert(0, str(E13 / "scripts")); import lib13 as L13  # noqa: E402
    lines = json.loads(subprocess.run(["git", "show", "HEAD:./analysis/lines.json"], cwd=E13, capture_output=True, check=True).stdout.decode("utf-8"))["lines"][::7]

    def amp(D, a, b):
        k = (TW >= a) & (TW <= b); live = np.abs(D).sum(0) > 0
        return np.sqrt(lib.smooth(np.mean(np.abs(np.fft.rfft(D[k][:, live] * np.hanning(k.sum())[:, None], NF, axis=0)) ** 2, axis=1), 5))
    r = np.array([np.log(amp(D, 2.3, 2.8) / amp(D, 1.2, 1.6)) for D in (L13.read_inline(ln) for ln in lines)]).mean(0); FIT = (F5 >= 10) & (F5 <= 55)
    res = lib.smooth(r - np.polyval(np.polyfit(F5[FIT], r[FIT], 1), F5), 11); i0, i1 = int(np.argmin(np.abs(F5 - 8))), int(np.argmin(np.abs(F5 - 72))); res[:i0] = res[i0]; res[i1:] = res[i1]
    CURVE.write_text(json.dumps(dict(f=F5.tolist(), r=res.tolist())), encoding="utf-8"); print("r(f) at 8..72 step 8:", " ".join("%+.2f" % res[int(np.argmin(np.abs(F5 - f)))] for f in range(8, 73, 8)))

if __name__ == "__main__" and sys.argv[1] == "run":
    cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); fp = float(cm["fpeak"]); out = ROOT / "results" / "probe2"; out.mkdir(parents=True, exist_ok=True)
    r = np.asarray(json.loads(CURVE.read_text(encoding="utf-8"))["r"]); F = np.fft.rfftfreq(1024, DT)
    for did in sys.argv[2:]:
        t0 = time.monotonic(); arr, meta = load(did); ideal, R_ = arr["ideal"], arr["R"]; rec = []
        for var, D in (("O", arr["D"]), ("F", tvfilt(arr["signal"], r) + arr["noise"])):
            w_out = sc11.prepare_sc(D, cm)["wv"]["w_out"]; inv = snr_rows(D, TW, F); rec.append(dict(var=var, method="Input", m=1.0, S=score(D, ideal)["S"], T2=t2(D, R_, fp), hf=hf(D)))
            for m in MS:
                o = CW.methods(D, cm, Q, m * meta["c_eff_cum"], TW, 1.4, w_out, inv); o["UQWS"] = lib12.tvw(o["UQW"], w_out, SEL["B3"]["eps"], win=SEL["B3"]["win"])
                rec += [dict(var=var, method=k, m=m, S=score(o[k], ideal)["S"], T2=t2(o[k], R_, fp), hf=hf(o[k])) for k in KEYS]
        (out / f"{did}.json").write_text(json.dumps(rec, indent=1), encoding="utf-8"); print(did, "%.0fs" % (time.monotonic() - t0), flush=True)
