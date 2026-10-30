"""Shared code of the SC-a ceiling check (PREDICTIONS_ZH.md): why SC-a (and the Wiener baseline VDB) stay below
70 dB inverse Q on the experiment-15 data even with the true effective Q.

Everything is imported unchanged from experiment 15 (lib15 and, through it, the experiment-7 to -13 code).
New code: a variant of lib.wavelets whose shallow band edge does not stop at a notch of the coherence SNR.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "fault_q_confirm_20261030" / "scripts"))
from lib15 import np, lib, semi, sc11, lib12, TW, DT, SEL, A_TUNED, G_IQ, META, load, score  # noqa: E402

IDS = [f"Z{k}_c{c}_r{r}" for k in "ABC" for c in (0.75, 1.05, 1.4) for r in (0, 1)]
E15 = BASE / "fault_q_confirm_20261030" / "results"


def band_last(f, snr, fpk):
    """Highest frequency at or above fpk where SNR >= 1 (notches allowed); lib.snr_band stops at the first notch."""
    ok = np.where((f >= f[int(np.argmin(np.abs(f - fpk)))]) & (snr >= 1))[0]
    return float(f[ok[-1]]) if len(ok) else float(f[int(np.argmin(np.abs(f - fpk)))])


def wavelets(D, cm, robust, nfft=512):
    """lib.wavelets; with robust=True the shallow edge fs uses band_last, everything else unchanged
    (the target edge fd keeps lib.snr_band)."""
    if not robust:
        return lib.wavelets(D, TW, cm, nfft)
    f, snr_s, P_s = lib.coherence_snr(D, TW, 1.20, 1.60, nfft)
    _, snr_d, _ = lib.coherence_snr(D, TW, 2.00, 3.00, nfft)
    amp = np.sqrt(P_s) / np.exp(-np.pi * f * lib.absorption(cm, [1.40])[0])
    fpk = float(f[np.argmax(np.sqrt(P_s) * (f >= 5))])
    fs = band_last(f, snr_s, fpk); fd = lib.snr_band(f, snr_d, fpk)

    def cut(fc):
        g = np.ones_like(f); r = (f > fc - 5) & (f < fc)
        g[r] = 0.5 * (1 + np.cos(np.pi * (f[r] - (fc - 5)) / 5)); g[f >= fc] = 0; g[f < 3] = 0
        return g

    def zp(a):
        w = np.roll(np.fft.irfft(a, nfft), 100)[:201]; return w / np.max(np.abs(w))

    a_in = amp * cut(fs); a_out = a_in * cut(min(fd, fs))
    w_in, w_out = zp(a_in), zp(a_out)
    w_out *= np.max(np.abs(np.fft.rfft(w_in, nfft))) / np.max(np.abs(np.fft.rfft(w_out, nfft)))
    return dict(w_in=w_in, w_out=w_out, lag=np.arange(-100, 101) * DT, f=f, snr_shallow=snr_s, snr_target=snr_d,
                f_peak=fpk, fmax_shallow=fs, fmax_target=fd)


def methods(D, cm, Q, c, p, wv0, wv1):
    """All with the true effective factor c = c_eff_cum. SC: a = 2^-13 (experiment 12 B2); mu_relative 10 or 0.
    0 = lib.wavelets as in experiment 15, 1 = robust shallow edge."""
    sc = lambda wv, mu: lib.solve(cm, Q * c, wv, D, p["sigma"], p["K"], p["live"], A_TUNED, mu)[0]
    vdb = lambda wv: lib12.tvw(lib12.iq(D, cm, Q, None, c, 0), wv["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"])
    out = dict(IQ70=lib12.iq(D, cm, Q, c, c, G_IQ), SC_0_mu10=sc(wv0, lib.SEL_SC[1]), VDB_0=vdb(wv0))
    out["SC_0_mu0"] = sc(wv0, 0.0)
    same = wv1["fmax_shallow"] == wv0["fmax_shallow"]
    out["SC_1_mu10"] = out["SC_0_mu10"] if same else sc(wv1, lib.SEL_SC[1])
    out["SC_1_mu0"] = out["SC_0_mu0"] if same else sc(wv1, 0.0)
    out["VDB_1"] = vdb(wv1)
    return out
