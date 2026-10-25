"""Freeze RESUME_ZH.md and the Octave-side code before any U1-U4 diagnosis runs."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
COMPAT = ROOT.parent / "octave_compat_20261025"
out = ROOT / "logs" / "RESUME_ADDENDUM_SHA256.json"
assert not out.exists() and not (ROOT / "predictions").exists() and not (ROOT / "cache" / "outputs").exists()
head = subprocess.run(["git", "-C", str(ROOT.parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
paths = [ROOT / "RESUME_ZH.md", *sorted((ROOT / "scripts_octave").glob("*")), *sorted(COMPAT.glob("*.m")), *sorted(COMPAT.glob("*.py")),
         *sorted((ROOT / "octave_check").glob("*.csv"))]
rec = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
           note="Before any confirmation diagnosis on U1-U4; predictions/ and cache/outputs/ do not exist.",
           hashes={p.relative_to(ROOT.parent).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()})
out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
print(f"{len(rec['hashes'])} files frozen at {head[:8]}")
