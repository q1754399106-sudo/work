"""Freeze PLAN_ZH.md, the scripts and the input snapshot before any feature is computed."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / "logs" / "PRE_RUN_HASHES.json"
assert not out.exists() and not list((ROOT / "analysis").glob("features_*.csv"))
head = subprocess.run(["git", "-C", str(ROOT.parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
compat = ROOT.parent / "octave_compat_20261025"
paths = [ROOT / "PLAN_ZH.md", ROOT / "cache" / "SNAPSHOT.json", *sorted((ROOT / "scripts").glob("*")),
         *sorted(compat.glob("*.m")), *sorted(compat.glob("*.py"))]
rec = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
           note="Before any feature computation.",
           hashes={p.relative_to(ROOT.parent).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()})
out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
print(f"{len(rec['hashes'])} files frozen")
