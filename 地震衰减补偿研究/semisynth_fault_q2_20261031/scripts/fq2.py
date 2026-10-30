"""Second batch of the fault / complex-Q semi-synthetic data (DATA_ZH.md): the generator of
semisynth_fault_q_20261030 (fq.py, unchanged) on three new structures D, E, F with new inlines and seeds.
Only the Z (complex Q) variant is kept; c0 = 0.75 and 1.05 (the attenuation closest to the real line)."""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; OLD = ROOT.parent / "semisynth_fault_q_20261030"
sys.path.insert(0, str(OLD / "scripts"))
import fq  # noqa: E402
from fq import np, semi  # noqa: E402,F401

STRUCT = {"D": 10253, "E": 10305, "F": 10351}
NOISE = {"D": (9991, 10259), "E": (9997, 10287), "F": (10331, 10365)}
GEO_SEED = {"D": 15004, "E": 15005, "F": 15006}
RC_SEED = {(k, r): 15130 + 10 * i + r for i, k in enumerate(STRUCT) for r in (0, 1)}
C0 = (0.75, 1.05)
fq.STRUCT.update(STRUCT); fq.NOISE.update(NOISE); fq.GEO_SEED.update(GEO_SEED); fq.RC_SEED.update(RC_SEED)
IDS = [f"Z{k}_c{c}_r{r}" for k in STRUCT for c in C0 for r in (0, 1)]


def load(did):
    """Arrays (float32, TW x NX) and metadata of a built dataset, e.g. load("ZD_c0.75_r0")."""
    meta = json.loads((ROOT / "analysis" / "datasets.json").read_text(encoding="utf-8"))[did]
    return dict(np.load(ROOT / "cache" / "data" / f"{did}.npz")), meta
