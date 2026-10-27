"""POST-HOC (not registered; RESULTS_ZH.md): is the boosted 50-90 Hz energy laterally coherent or not?

Gate A1 (60-100 Hz energy ratio <= 10) was mis-specified: correct compensation of Q ~ 100 at 2.5 s must raise
that ratio by thousands. This script reports, in the target window 2.0-3.0 s, the lateral-coherence SNR of
the input and of each output per frequency band, and the coherent / incoherent energy gain per band:
  coherent gain   = C_out / C_in      (C: mean real cross-power of adjacent traces)
  incoherent gain = (P_out - C_out) / (P_in - C_in)
Caveat: SC and the coherence Standard contain lateral regularisation, which by itself raises coherence of
the output; the ratios describe the outputs, they do not prove that the boosted energy is reflection signal.
Writes analysis/posthoc_coherence.json.
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib  # noqa: E402
from lib import np, bm  # noqa: E402
from scipy.io import loadmat  # noqa: E402
import run  # noqa: E402

BANDS = [(10, 30), (30, 50), (50, 70), (70, 90), (90, 110)]


def parts(A, tw):
    w = (tw >= 2.0 - 1e-9) & (tw <= 3.0 + 1e-9)
    X = np.fft.rfft(A[w] * lib.taper(w.sum())[:, None], 512, axis=0)
    P = np.mean(np.abs(X) ** 2, axis=1)
    C = np.mean(np.real(X[:, :-1] * np.conj(X[:, 1:])), axis=1)
    return np.fft.rfftfreq(512, lib.DT), P, C


def table(outs, tw):
    f, Pi, Ci = parts(outs["Input"], tw)
    res = {}
    for m, A in outs.items():
        _, P, C = parts(np.asarray(A, float), tw)
        rows = []
        for lo, hi in BANDS:
            b = (f >= lo) & (f < hi)
            snr = max(C[b].sum(), 0) / max(P[b].sum() - max(C[b].sum(), 0), 1e-300)
            rows.append(dict(band_Hz=[lo, hi], snr=float(snr),
                             coherent_gain=float(C[b].sum() / Ci[b].sum()),
                             incoherent_gain=float((P[b].sum() - C[b].sum()) / (Pi[b].sum() - Ci[b].sum()))))
        res[m] = rows
    return res


out = {}
for well, p in run.WELLS.items():
    c = loadmat(p["curtain"], simplify_cells=True); tw = bm.vec(c["tw"])
    outs = {"Input": np.asarray(c["D"], float), "Standard": loadmat(p["standard"], simplify_cells=True)["G"]["Standard"],
            "SC": loadmat(p["sc"], simplify_cells=True)["X"],
            "SCa": np.load(run.OUT / f"{well}_SCa.npz")["out"], "QSCa": np.load(run.OUT / f"{well}_QSCa.npz")["out"]}
    out[well] = table(outs, tw)
d = np.load(run.OUT / "IL10141_input.npz")
out["IL10141"] = table({"Input": d["D"], "SCa": np.load(run.OUT / "IL10141_SCa.npz")["out"],
                        "QSCa": np.load(run.OUT / "IL10141_QSCa.npz")["out"]}, d["tw"])
bm.writejson(run.AN / "posthoc_coherence.json", out)
for k, v in out.items():
    print(k)
    for m, rows in v.items():
        print(f"  {m:9}", "  ".join(f"{r['band_Hz'][0]}-{r['band_Hz'][1]}: SNR {r['snr']:5.2f} coh x{r['coherent_gain']:8.1f} "
                                      f"inc x{r['incoherent_gain']:8.1f}" for r in rows[1:4]))
