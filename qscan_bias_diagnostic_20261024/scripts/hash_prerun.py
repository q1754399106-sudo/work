"""Record SHA-256 of the frozen plan, manifest, scripts and input snapshot before the main run."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
head = subprocess.run(["git", "-C", str(ROOT.parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
paths = [ROOT / "PLAN_ZH.md", ROOT / "predictions" / "manifest.csv", ROOT / "cache" / "SNAPSHOT.json",
         *sorted((ROOT / "scripts").glob("*.m")), *sorted((ROOT / "scripts").glob("*.py"))]
out = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
           note="Before diagnosing any manifest gather; development probe, input checks, smoke and timing tests only.",
           hashes={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(ROOT / "logs" / "PRE_RUN_HASHES.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(len(paths), "files hashed")
