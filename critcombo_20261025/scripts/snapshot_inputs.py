"""Copy every Git-archived input this round reads into the ignored cache/ and verify it byte for
byte, and record SHA-256 of the sweep-cache replicates (r001-r006) the scoring reads.

Sources are Git objects of the current HEAD (regular blobs, or local Git LFS objects checked
against their pointer oid), so the run does not depend on the shared sparse checkout, which
other sessions rewrite. Writes cache/SNAPSHOT.json. (qscan_offset_20261025/scripts/snapshot_inputs.py
with this round's file list and the sweep-cache hashes.)
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DEST = ROOT / "cache" / "inputs"


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=True).stdout


def blob(path):
    data = git("show", f"HEAD:{path}")
    if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
        fields = dict(line.split(" ", 1) for line in data.decode().strip().splitlines())
        oid, size = fields["oid"].split(":", 1)[1], int(fields["size"])
        data = (REPO / ".git" / "lfs" / "objects" / oid[:2] / oid[2:4] / oid).read_bytes()
        assert len(data) == size and hashlib.sha256(data).hexdigest() == oid, path
    return data


head = git("rev-parse", "HEAD").decode().strip()
paths = git("ls-tree", "-r", "--name-only", "HEAD", "bae_q_20261024/core", "bae_q_20261024/source",
            "bae_q_20261024/inputs/independent", "bae_q_20261024/truth/independent",
            "qscan_offset_20261025/models", "qscan_offset_20261025/truth").decode().split()
files = {}
for path in paths:
    data = blob(path)
    out = DEST / path
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    files[path] = dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
sweeps = {}
for base in ("combo_bias_variance_20261024", "qscan_offset_20261025"):
    for f in sorted((REPO / base / "cache" / "sweep").glob("*/r00[1-6].mat")):
        sweeps[f.relative_to(REPO).as_posix()] = hashlib.sha256(f.read_bytes()).hexdigest()
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                             commit=head, files=files, sweep_cache=sweeps), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files from {head[:8]}, {len(sweeps)} sweep replicates hashed")
