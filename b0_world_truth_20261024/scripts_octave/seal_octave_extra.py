"""Supplementary seal (RESUME_ZH.md): Octave-side files that the unchanged scripts/seal_confirm.py
does not cover. Run right after seal_confirm.py, before any truth is read."""
from pathlib import Path
import datetime, hashlib, json

ROOT = Path(__file__).resolve().parents[1]
COMPAT = ROOT.parent / "octave_compat_20261025"
P = ROOT / "predictions"
assert (P / "SEALED_SHA256.json").exists() and not (P / "SEALED_OCTAVE_SHA256.json").exists()
assert not list((ROOT / "analysis").glob("truth_U*.csv")) and not list((ROOT / "analysis").glob("combo_truth_U*.csv"))
paths = [ROOT / "RESUME_ZH.md", *sorted((ROOT / "scripts_octave").glob("*")), *sorted(COMPAT.glob("*.m")), *sorted(COMPAT.glob("*.py"))]
receipt = dict(sealed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), truth_read_before_seal=False,
               hashes={p.relative_to(ROOT.parent).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()})
(P / "SEALED_OCTAVE_SHA256.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
print(f"sealed {len(receipt['hashes'])} Octave-side files")
