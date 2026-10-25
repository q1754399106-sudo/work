"""Record SHA-256 of the plan, every script and the input snapshot (which holds the sweep-cache hashes) before any criterion runs."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / "logs" / "PRE_RUN_HASHES.json"
assert not out.exists(), "Pre-run hashes exist; refusing to overwrite."
assert not (ROOT / "cache" / "crit").exists(), "Criterion outputs exist; the plan must be frozen first."
paths = ["PLAN_ZH.md", "cache/SNAPSHOT.json"] + sorted(
    p.relative_to(ROOT).as_posix() for p in (ROOT / "scripts").iterdir() if p.suffix in {".m", ".py"})
head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
out.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
                               note="Before any criterion computation.", hashes=hashes), indent=2), encoding="utf-8")
print(f"{len(hashes)} hashes written")
