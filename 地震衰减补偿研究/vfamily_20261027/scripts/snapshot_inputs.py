"""Build cache/runner (round3 core/, source/ and run_combo.m from HEAD, byte for byte) and
cache/runner/run_combo3.m, which differs from run_combo.m only in the function name and the
Q-error list {-30, 0, +30}; each replaced text must occur exactly once. Copies the layered-scene
config (soft_refit primary.json) and the calibration policy (round3 candidate_policy.json) used for
V4/V5, and the novelty SC code and parameters at ff95a69 (checked equal to HEAD). Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
REPO = RESEARCH.parent
RUNNER = ROOT / "cache" / "runner"
NOVELTY = ["novelty_evidence_20261025/scripts/benchmark.py", "novelty_evidence_20261025/scripts/common.py",
           "novelty_evidence_20261025/analysis/selection.json"]
SPECS = ["adaptive_cv_plan.json", "candidate_policy.json", "options.json", "p1_experiment_spec.json",
         "stability_spec.json", "uq_spec.json", "v2_options.json"]
REPLACE = [("function run_combo(sets)", "function run_combo3(sets)"),
           ("for e=[-50 -30 -10 0 10 30 50]", "for e=[-30 0 30]")]


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=True).stdout


assert not (RUNNER / "dev").exists() and not (RUNNER / "results").exists(), "runner already used"
files = {}
r3 = "地震衰减补偿研究/round3_20261024"
for path in git("ls-tree", "-r", "--name-only", "HEAD", f"{r3}/core", f"{r3}/source", f"{r3}/run_combo.m").decode().split():
    data = git("show", f"HEAD:{path}")
    assert not data.startswith(b"version https://git-lfs"), path
    out = RUNNER / path[len(r3) + 1:]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    files[path] = dict(commit="HEAD", sha256=hashlib.sha256(data).hexdigest())
# Specs the calibration chain reads from <root>/inputs (same blobs as coherence_soft_20261023/inputs,
# the root V5 was calibrated from).
for b in SPECS:
    path = f"{r3}/inputs/{b}"
    data = git("show", f"HEAD:{path}")
    assert data == git("show", f"HEAD:地震衰减补偿研究/coherence_soft_20261023/inputs/{b}"), b
    out = RUNNER / "inputs" / b
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    files[path] = dict(commit="HEAD", sha256=hashlib.sha256(data).hexdigest())
text = (RUNNER / "run_combo.m").read_bytes().decode("utf-8")
for old, new in REPLACE:
    assert text.count(old) == 1, old
    text = text.replace(old, new)
(RUNNER / "run_combo3.m").write_bytes(text.encode("utf-8"))
files["cache/runner/run_combo3.m"] = dict(derived_from=f"{r3}/run_combo.m", replacements=REPLACE,
                                           sha256=hashlib.sha256(text.encode("utf-8")).hexdigest())
for path in ("地震衰减补偿研究/soft_refit_20261023/inputs/primary.json", f"{r3}/inputs/candidate_policy.json"):
    data = git("show", f"HEAD:{path}")
    out = ROOT / "cache" / "inputs" / Path(path).name
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    files[path] = dict(commit="HEAD", sha256=hashlib.sha256(data).hexdigest())
for path in NOVELTY:
    data = git("show", f"ff95a69:{path}")
    assert data == git("show", f"HEAD:地震衰减补偿研究/{path}"), f"{path} changed after ff95a69"
    out = ROOT / "cache" / "inputs" / path
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    files[path] = dict(commit="ff95a69", sha256=hashlib.sha256(data).hexdigest())
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), files=files), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files")
