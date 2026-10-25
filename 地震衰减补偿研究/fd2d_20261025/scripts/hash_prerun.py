"""Freeze points of PLAN_ZH.md.

hash_prerun.py plan    after the verification, before any H1-H6 shot: plan, scripts, input
                       snapshot and verification results -> logs/PRE_RUN_HASHES.json
hash_prerun.py models  after assemble_models, before any sweep replicate: checks that every
                       file frozen at 'plan' is unchanged and adds models/ and truth/
                       -> logs/PRE_SWEEP_HASHES.json
"""
from pathlib import Path
import datetime, hashlib, json, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]
mode = sys.argv[1]
assert mode in {"plan", "models"}
sha = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
frozen = ["PLAN_ZH.md", "cache/SNAPSHOT.json", "analysis/verification.json"] + sorted(
    p.relative_to(ROOT).as_posix() for p in (ROOT / "scripts").iterdir() if p.suffix in {".m", ".py"})
if mode == "plan":
    out = ROOT / "logs" / "PRE_RUN_HASHES.json"
    assert not out.exists(), "Plan hashes exist; refusing to overwrite."
    assert not (ROOT / "cache" / "shots").exists() and not (ROOT / "models").exists(), "Shots or models exist already."
    hashes = {p: sha(p) for p in frozen}
    note = "After the forward verification (test geologies only), before any H1-H6 shot."
else:
    out = ROOT / "logs" / "PRE_SWEEP_HASHES.json"
    assert not out.exists(), "Sweep hashes exist; refusing to overwrite."
    assert not (ROOT / "cache" / "sweep").exists(), "Sweep outputs exist; models must be frozen first."
    plan = json.loads((ROOT / "logs" / "PRE_RUN_HASHES.json").read_text(encoding="utf-8"))["hashes"]
    changed = [p for p, h in plan.items() if sha(p) != h]
    assert not changed, f"Frozen files changed: {changed}"
    hashes = {p: sha(p) for p in sorted(
        q.relative_to(ROOT).as_posix() for d in ("models", "truth") for q in (ROOT / d).iterdir() if q.suffix == ".mat")}
    note = "After assemble_models (H1-H6 generated), before any sweep replicate; plan-time hashes re-verified."
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), git_head=head,
                               note=note, hashes=hashes), indent=2), encoding="utf-8")
print(f"{mode}: {len(hashes)} hashes written to {out.name}")
