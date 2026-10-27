"""Restore every input named in PLAN_ZH.md into the ignored cache/inputs and verify it.

Git files come from their registered commits (LFS resolved through `git lfs smudge`, checked against the
pointer oid) and must equal the HEAD blob. The MIG_RAW cube is the local copy restored for experiment 6,
checked against its SHA-256. Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
PREFIX = "地震衰减补偿研究/"
DEST = ROOT / "cache" / "inputs"
E6 = "field_vsp_volve_f11_20261026/"
E5 = "field_vsp_volve_20261025/"
SOURCES = {  # commit -> (paths at that commit, prefixed?)
    "ff95a69": (["novelty_evidence_20261025/scripts/benchmark.py", "novelty_evidence_20261025/scripts/common.py",
                 "novelty_evidence_20261025/analysis/selection.json"], False),
    "609ddcd": ([E5 + "inputs/curtain.mat", E5 + "analysis/references.mat", E5 + "analysis/vsp_q.json",
                 E5 + "results/V0.mat"], False),
    "0b04ea9": ([E6 + "inputs/posthoc_curtain.mat", E6 + "analysis/posthoc_references.mat",
                 E6 + "analysis/posthoc_vsp_q.json", E6 + "results_posthoc/V0.mat",
                 E6 + "results_sc_posthoc/F11_V0_SC.mat", E6 + "results_sc_posthoc/F11_V0_QSC.mat",
                 E6 + "results_sc/F15A_V0_SC.mat", E6 + "results_sc/F15A_V0_QSC.mat"], True),
}
CUBE = (ROOT.parent / E6 / "cache" / "seismic" / "MIG_RAW.segy",
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
for rev, (paths, prefixed) in SOURCES.items():
    for path in paths:
        data = blob(rev, PREFIX + path if prefixed else path)
        assert data == blob("HEAD", PREFIX + path), f"{path} changed after {rev}"
        out = DEST / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[path] = dict(commit=rev, size=len(data), sha256=hashlib.sha256(data).hexdigest())
assert hashlib.sha256(CUBE[0].read_bytes()).hexdigest() == CUBE[1], "MIG_RAW cube differs"
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), cube=str(CUBE[0]), cube_sha256=CUBE[1], files=files),
    indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files")
