"""SC-a at 4 x noise with the notch-tolerant wavelet band rule (exploratory, NOTE_ZH.md section 11.3).

At 4 x noise lib.wavelets stops the shallow band edge at the first notch of the coherence SNR (36-61 Hz instead of
98-125 Hz at calibrated noise): the mechanism 1 of sca_ceiling_20261030, whose fix (band_last for the shallow edge,
libc.wavelets robust=True) is copied here unchanged. Variant 'both' also applies band_last to the target edge
(descriptive; the repo never used it). SR route, the 16 evaluation datasets of pilot_sc_tune.py.

  R1_sel   robust shallow edge, (a, mu) selected by pilot_sc_tune.py (2^-9, 2)
  R1_def   robust shallow edge, a 2^-13, mu 0.25 (experiments 12 / 17)
  R2_sel   both edges notch-tolerant, (2^-9, 2)

usage: pilot_sc_robust.py run K N
"""
import sys, json, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import uqp  # noqa: E402,F401
from uqp import np, lib, semi, ROOT, DT  # noqa: E402
from uqw import load, TW  # noqa: E402
from lib17 import score, sc11  # noqa: E402
from pilot_sc_tune import EVAL  # noqa: E402

OUT = ROOT / "results" / "pilot_sc_robust"
VARIANTS = {"R1_sel": (False, 2.0 ** -9, 2.0), "R1_def": (False, 2.0 ** -13, 0.25), "R2_sel": (True, 2.0 ** -9, 2.0)}


def band_last(f, snr, fpk):
    """sca_ceiling_20261030/scripts/libc.py, unchanged."""
    ok = np.where((f >= f[int(np.argmin(np.abs(f - fpk)))]) & (snr >= 1))[0]
    return float(f[ok[-1]]) if len(ok) else float(f[int(np.argmin(np.abs(f - fpk)))])


def wavelets(D, cm, both, nfft=512):
    """libc.wavelets(robust=True); with both=True the target edge fd also uses band_last."""
    f, snr_s, P_s = lib.coherence_snr(D, TW, 1.20, 1.60, nfft)
    _, snr_d, _ = lib.coherence_snr(D, TW, 2.00, 3.00, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * lib.absorption(cm, [1.40])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))])
    fs = band_last(f, snr_s, fpk); fd = band_last(f, snr_d, fpk) if both else lib.snr_band(f, snr_d, fpk)

    def cut(fc):
        g = np.ones_like(f); r = (f > fc - 5) & (f < fc)
        g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fc - 5)) / 5)); g[f >= fc] = 0; g[f < 3] = 0
        return g

    def zp(a):
        w = np.roll(np.fft.irfft(a, nfft), 100)[:201]; return w / np.max(np.abs(w))

    a_in = amp * cut(fs); a_out = a_in * cut(min(fd, fs)); w_in, w_out = zp(a_in), zp(a_out)
    w_out *= np.max(np.abs(np.fft.rfft(w_in, nfft))) / np.max(np.abs(np.fft.rfft(w_out, nfft)))
    return dict(w_in=w_in, w_out=w_out, fmax_shallow=fs, fmax_target=fd)


if __name__ == "__main__":
    k, n = int(sys.argv[2]), int(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True); cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)
    jobs = [(d, v) for v in VARIANTS for d in EVAL]; cache = {}
    for j, (did, v) in enumerate(jobs):
        path = OUT / f"{did}_{v}.json"
        if j % n != k or path.exists():
            continue
        both, a, mu = VARIANTS[v]; t0 = time.monotonic()
        if did not in cache:
            cache.clear(); arr, meta = load(did); D = arr["signal"] + 4.0 * arr["noise"]
            c = json.loads((ROOT / "results" / "pilot_semi" / f"{did}_nd.json").read_text(encoding="utf-8"))["sr"]["c_sr"]
            cache[did] = (arr, D, c, sc11.prepare_sc(D, cm))
        arr, D, c, p = cache[did]; wv = wavelets(D, cm, both)
        X, info, _ = lib.solve(cm, Q * c, wv, D, p["sigma"], p["K"], p["live"], a, mu)
        out = dict(id=did, variant=v, a=a, mu=mu, c_sr=c, fmax_shallow=wv["fmax_shallow"], fmax_target=wv["fmax_target"],
                   S=float(score(X, arr["ideal"])["S"]), converged=bool(info["converged"]), seconds=time.monotonic() - t0)
        path.write_text(json.dumps(out, indent=1), encoding="utf-8"); print(out, flush=True)
