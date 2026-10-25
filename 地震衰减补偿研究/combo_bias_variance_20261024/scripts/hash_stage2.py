"""Record SHA-256 of the stage-2 plan and every script before the sweep runs (only the smoke replicate may exist)."""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / "logs" / "PRE_STAGE2_HASHES.json"
assert not out.exists(), "Stage-2 hashes exist; refusing to overwrite."
done = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "cache" / "sweep").rglob("r*.mat"))
assert done == ["cache/sweep/S01_F1_Q100_14dB/r001.mat"], done
paths = ["PLAN_STAGE2_ZH.md", "PLAN_ZH.md", "cache/SNAPSHOT.json"] + sorted(
    p.relative_to(ROOT).as_posix() for p in (ROOT / "scripts").iterdir() if p.suffix in {".m", ".py"})
head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
out.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
                               note="After the one-replicate smoke check (S01 r=1, reproduction only), before the sweep.",
                               hashes=hashes), indent=2), encoding="utf-8")
print(f"{len(hashes)} hashes written")
