"""Freeze PLAN_ZH.md, scripts, compatibility layer and input snapshot before K1-K8 exist."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / "logs" / "PRE_RUN_HASHES.json"
assert not out.exists() and not (ROOT / "models").exists() and not (ROOT / "truth").exists()
head = subprocess.run(["git", "-C", str(ROOT.parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
compat = ROOT.parent / "octave_compat_v2_20261025"
paths = [ROOT / "PLAN_ZH.md", ROOT / "cache" / "SNAPSHOT.json", ROOT / "analysis" / "generator_check.json",
         *sorted((ROOT / "scripts").glob("*")), *sorted(compat.glob("*.m")), *sorted(compat.glob("*.py"))]
rec = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
           note="Before generating K1-K8; no confirmation data exist.",
           hashes={p.relative_to(ROOT.parent).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
print(f"{len(paths)} files frozen at {head[:8]}")
