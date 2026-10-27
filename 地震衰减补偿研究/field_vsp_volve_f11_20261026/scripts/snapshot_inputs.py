"""Copy every frozen input named in PLAN_ZH.md into the ignored cache/ and verify it.

Code and experiment-5 files are read from their registered commits (those commits predate the move into
地震衰减补偿研究/, so the path there has no prefix) and must equal the HEAD blob under the prefix.
LFS files are resolved through `git lfs smudge` and checked against the pointer oid.
The MIG_RAW cube restored into cache/seismic/ is checked against its SHA-256.
Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
PREFIX = "地震衰减补偿研究/"
DEST = ROOT / "cache" / "inputs"
SOURCES = {  # registered commit -> paths (PLAN_ZH.md, "方法")
    "b59a6c3": ["bae_q_20261024/core", "bae_q_20261024/source"],
    "f38f3e3": ["round3_20261024/core/coherence_scores.m", "round3_20261024/core/apply_gain.m"],
    "16aabfe": ["combo_bias_variance_20261024/scripts/combo_outputs.m",
                "combo_bias_variance_20261024/scripts/qscan_select.m"],
    "ff95a69": ["novelty_evidence_20261025/scripts/benchmark.py", "novelty_evidence_20261025/scripts/common.py",
                "novelty_evidence_20261025/analysis/selection.json"],
    "609ddcd": ["field_vsp_volve_20261025/inputs/curtain.mat", "field_vsp_volve_20261025/analysis/vsp_q.json",
                "field_vsp_volve_20261025/analysis/references.mat", "field_vsp_volve_20261025/analysis/scores.csv",
                "field_vsp_volve_20261025/results"],
}
CUBE = (ROOT / "cache" / "seismic" / "MIG_RAW.segy",
        "22dbd614ee34a8d2becf3ea66d55db6887b8ed7804212452077aef9340c9ef92")


def git(*args, stdin=None):
    return subprocess.run(["git", "-C", str(REPO), *args], input=stdin, capture_output=True, check=True).stdout


def blob(rev, path):
    data = git("show", f"{rev}:{path}")
    if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
        oid = dict(line.split(" ", 1) for line in data.decode().strip().splitlines())["oid"].split(":", 1)[1]
        data = git("lfs", "smudge", "--", path, stdin=data)
        assert hashlib.sha256(data).hexdigest() == oid, path
    return data


files = {}
for rev, paths in SOURCES.items():
    for path in git("ls-tree", "-r", "--name-only", rev, *paths).decode().split():
        data = blob(rev, path)
        assert data == blob("HEAD", PREFIX + path), f"{path} changed after {rev}"
        out = DEST / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[path] = dict(commit=rev, size=len(data), sha256=hashlib.sha256(data).hexdigest())
assert hashlib.sha256(CUBE[0].read_bytes()).hexdigest() == CUBE[1], "MIG_RAW cube differs from experiment 5"
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), cube_sha256=CUBE[1], files=files), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files")
