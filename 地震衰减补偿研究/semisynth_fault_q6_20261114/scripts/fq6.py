"""Sixth batch of the fault / complex-Q semi-synthetic data (DATA_ZH.md): the generator of
semisynth_fault_q_20261030 (fq.py, unchanged) on three new structures P, R, T (letters Q and S skipped: they name the quality factor and the score) with new inlines and seeds.
Only the Z (complex Q) variant is kept; c0 = 0.75 and 1.05 (the attenuation closest to the real line)."""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; OLD = ROOT.parent / "semisynth_fault_q_20261030"
sys.path.insert(0, str(OLD / "scripts"))
import fq  # noqa: E402
from fq import np, semi  # noqa: E402,F401

STRUCT = {"P": 10023, "R": 10089, "T": 10197}
NOISE = {"P": (10009, 10055), "R": (10079, 10131), "T": (10167, 10211)}
GEO_SEED = {"P": 20261141, "R": 20261142, "T": 20261143}
RC_SEED = {(k, r): 20261150 + 10 * i + r for i, k in enumerate(STRUCT) for r in (0, 1)}
C0 = (0.75, 1.05)
fq.STRUCT.update(STRUCT); fq.NOISE.update(NOISE); fq.GEO_SEED.update(GEO_SEED); fq.RC_SEED.update(RC_SEED)
IDS = [f"Z{k}_c{c}_r{r}" for k in STRUCT for c in C0 for r in (0, 1)]


def load(did):
    """Arrays (float32, TW x NX) and metadata of a built dataset, e.g. load("ZP_c0.75_r0")."""
    meta = json.loads((ROOT / "analysis" / "datasets.json").read_text(encoding="utf-8"))[did]
    return dict(np.load(ROOT / "cache" / "data" / f"{did}.npz")), meta
