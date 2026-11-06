"""Post hoc look at T4 from stored results only: which model and lambda drive the waveform of fc* under the discrepancy lambda?"""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"
D = json.loads((OUT / "test_disc.json").read_text(encoding="utf-8")); R = json.loads((OUT / "test.json").read_text(encoding="utf-8"))
for qf in (0.8, 1.0, 1.2):
    for M in ("M4", "M5"):
        print(f"Q x{qf} {M}: " + "; ".join(f"{m} lam {next(r['lam_med'] for r in D if r['M'] == M and r['qf'] == qf and r['m'] == m):.2e} wave {next(r['wave'] for r in D if r['M'] == M and r['qf'] == qf and r['m'] == m):.3f} |bias| {np.mean([abs(r['meas'] - r['truth']) for r in D if r['M'] == M and r['qf'] == qf and r['m'] == m]):.3f}" for m in ("fc_star", "split", "direct")))
print("\nFixed lambda, 18 dB, Q correct: waveform by model")
for M in ("M4", "M5"):
    for lam in (1e-5, 1e-4, 1e-3, 1e-2):
        print(f"{M} {lam:.0e}: " + "; ".join(f"{m} {next(r['wave'] for r in R if r['M'] == M and r['snr'] == 18 and r['qf'] == 1.0 and r['lam'] == lam and r['m'] == m):.3f}" for m in ("fc_star", "split", "direct")))
