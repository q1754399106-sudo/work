"""Sixth batch of the fault / complex-Q semi-synthetic data (DATA_ZH.md; uqw_shaping_confirm_20261119/PLAN_CONFIRM_SHAPING_ZH.md
section 3.2-3.3): the generator of semisynth_fault_q_20261030 (fq.py, unchanged) on six new structures P R S T U V with the
lines of uqw_shaping_confirm_20261119/lines6.json, Z variant only, c0 = 0.75 and 1.05, r = 0 and 1.

New here: a non-constant-Q spectral departure b(f) added to the signal only (noise, ideal, R unchanged):
  h(f) = -min(max(f, 8), f_k); b = kappa [h - straight-line fit of h over 10-55 Hz]; end values below 8 / above 72 Hz;
  11-bin smoothing on the NF = 512 grid; kappa such that the 10-55 Hz curvature of b (specform.py: quadratic coefficient
  x (45/2)^2) is K_B = +0.57; f_k = default_rng(geology seed + 500).uniform(30, 42) per structure.
  TV: probe2.tvfilt(signal, b) (log gain b (t - 1.4) / 1.15);  ST: time-invariant gain exp(b) on the whole trace.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
sys.path.insert(0, str(BASE / "semisynth_fault_q_20261030" / "scripts"))
import fq  # noqa: E402
from fq import np, semi, lib, TW, DT  # noqa: E402,F401

STRUCT = {"P": 10023, "R": 10047, "S": 10073, "T": 10123, "U": 10167, "V": 10191}
NOISE = {"P": (10009, 10035), "R": (10055, 10063), "S": (10077, 10093), "T": (10117, 10143), "U": (10155, 10163), "V": (10179, 10211)}
GEO_SEED = {k: 15016 + i for i, k in enumerate(STRUCT)}
RC_SEED = {(k, r): 15250 + 10 * i + r for i, k in enumerate(STRUCT) for r in (0, 1)}
C0 = (0.75, 1.05); RS = (0, 1); K_B = 0.57; NF = 512; F5 = np.fft.rfftfreq(NF, DT); FIT = (F5 >= 10) & (F5 <= 55)
fq.STRUCT.update(STRUCT); fq.NOISE.update(NOISE); fq.GEO_SEED.update(GEO_SEED); fq.RC_SEED.update(RC_SEED)
IDS = [f"Z{k}_c{c}_r{r}" for k in STRUCT for c in C0 for r in RS]


def curvature(y, f=F5, fit=FIT):
    """10-55 Hz curvature as specform.py: quadratic coefficient x (45/2)^2."""
    return float(np.polyfit(f[fit], y[fit], 2)[0] * (45 / 2) ** 2)


def f_knee(key):
    return float(np.random.default_rng(GEO_SEED[key] + 500).uniform(30, 42))


def bend(key):
    """b(f) on F5 for structure key (kappa scaled to curvature K_B)."""
    h = -np.minimum(np.maximum(F5, 8.0), f_knee(key)); u = h - np.polyval(np.polyfit(F5[FIT], h[FIT], 1), F5)
    i0, i1 = int(np.argmin(np.abs(F5 - 8))), int(np.argmin(np.abs(F5 - 72))); u[:i0] = u[i0]; u[i1:] = u[i1]
    u = lib.smooth(u, 11); return u * (K_B / curvature(u))


def st_filter(X, b, nfft=1024):
    """Time-invariant zero-phase gain exp(b(f)) on whole traces (b interpolated to the nfft grid)."""
    n = X.shape[0]; g = np.exp(np.interp(np.fft.rfftfreq(nfft, DT), F5, b))
    return np.fft.irfft(np.fft.rfft(X, nfft, axis=0) * g[:, None], nfft, axis=0)[:n]


def load(did):
    """Arrays (float32, TW x NX: D, signal, signal_TV, signal_ST, noise, ideal, src, R, tstar, invQ, ...) and metadata."""
    meta = json.loads((ROOT / "analysis" / "datasets.json").read_text(encoding="utf-8"))[did]
    return dict(np.load(ROOT / "cache" / "data" / f"{did}.npz")), meta
