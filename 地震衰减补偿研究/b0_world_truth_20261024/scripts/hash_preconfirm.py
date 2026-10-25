"""Record SHA-256 of the frozen confirmation plan, development result and all code before U1-U4 exist."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
assert not (ROOT / "inputs").exists() and not (ROOT / "truth").exists(), "Unseen models already exist."
head = subprocess.run(["git", "-C", str(ROOT.parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
paths = [ROOT / "PLAN_ZH.md", ROOT / "DEV_PLAN_ZH.md", ROOT / "analysis" / "dev_selection.json", ROOT / "cache" / "SNAPSHOT.json",
         *sorted((ROOT / "analysis").glob("dev_w*.csv")), *sorted((ROOT / "scripts").glob("*.m")), *sorted((ROOT / "scripts").glob("*.py"))]
out = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
           note="Before generating U1-U4; no confirmation data existed.",
           hashes={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(ROOT / "logs" / "PRE_CONFIRM_HASHES.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(len(paths), "files hashed")
