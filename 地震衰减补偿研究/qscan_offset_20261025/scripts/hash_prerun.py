"""Record SHA-256 of the plan, every script, the input snapshot and the generated models before any sweep runs."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / "logs" / "PRE_RUN_HASHES.json"
assert not out.exists(), "Pre-run hashes exist; refusing to overwrite."
assert not (ROOT / "cache" / "sweep").exists(), "Sweep outputs exist; the plan must be frozen first."
paths = ["PLAN_ZH.md", "cache/SNAPSHOT.json"] + sorted(
    p.relative_to(ROOT).as_posix() for d in ("scripts", "models", "truth") for p in (ROOT / d).iterdir() if p.suffix in {".m", ".py", ".mat"})
head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
out.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
                               note="After make_dev (G1-G8 generated), before any sweep replicate.", hashes=hashes), indent=2), encoding="utf-8")
print(f"{len(hashes)} hashes written")
