"""Post-hoc diagnostic (not registered; written after seeing that C1 is degenerate): does the
shallow/deep consistency premise of the Cheng criterion hold for the elastic targets themselves?
For each model's target T (the scoring target, same domain as C = W*(X/scale)) the window spectra
are computed exactly as cheng_stats.m (thirds and halves, demean, Hann, nfft 512, 5-80 Hz), and
reported as the deep/shallow band-RMS amplitude ratio and the C-A / C-C relative error of T.
Writes analysis/posthoc_truth.csv.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "novelty_evidence_20261025" / "cache" / "deps"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.io import loadmat  # noqa: E402

IN = ROOT / "cache" / "inputs"
NFFT, DT = 512, 0.004
WIN = {"thirds": (slice(0, 150), slice(301, 451)), "halves": (slice(0, 225), slice(226, 451))}
f = np.arange(NFFT // 2 + 1) / (NFFT * DT)
band = (f >= 5) & (f <= 80)


def spectrum(C):
    L = C.shape[0]
    h = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(L) / (L - 1))
    F = np.abs(np.fft.fft((C - C.mean(0)) * h[:, None], NFFT, axis=0))[: NFFT // 2 + 1]
    return F.mean(1)[band]


rows = []
for d, n, g in (("qscan_offset_20261025", 8, "G"), ("fd2d_20261025", 6, "H")):
    for m in range(1, n + 1):
        T = loadmat(IN / d / "truth" / f"{g}{m}_truth.mat")["target"]
        row = dict(model=f"{g}{m}")
        for w, (s, dp) in WIN.items():
            As, Ad = spectrum(T[s]), spectrum(T[dp])
            row[f"{w}_deep_over_shallow_rms"] = float(np.linalg.norm(Ad) / np.linalg.norm(As))
            row[f"{w}_relative_error"] = float(np.linalg.norm(Ad - As) / np.linalg.norm(As))
        rows.append(row)
out = pd.DataFrame(rows)
out.to_csv(ROOT / "analysis" / "posthoc_truth.csv", index=False, float_format="%.4f")
print(out.round(3).to_string(index=False))
