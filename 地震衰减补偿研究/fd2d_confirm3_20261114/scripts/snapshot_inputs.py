"""Copy every input of the FD forward modelling into the ignored cache/ and verify it byte for byte.

fd2d_20261025/scripts/snapshot_inputs.py with the same file list; only the Git paths carry the research-folder prefix
(the folder moved under 地震衰减补偿研究/ on 2026-10-26) and LFS objects are read from the repository's .git/lfs.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]; REPO = ROOT.parents[1]; PREFIX = "地震衰减补偿研究/"; DEST = ROOT / "cache" / "inputs"


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), "-c", "core.quotepath=off", *args], capture_output=True, check=True).stdout


def blob(path):
    data = git("show", f"HEAD:{PREFIX}{path}")
    if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
        fields = dict(line.split(" ", 1) for line in data.decode().strip().splitlines())
        oid, size = fields["oid"].split(":", 1)[1], int(fields["size"])
        data = (REPO / ".git" / "lfs" / "objects" / oid[:2] / oid[2:4] / oid).read_bytes()
        assert len(data) == size and hashlib.sha256(data).hexdigest() == oid, path
    return data


head = git("rev-parse", "HEAD").decode().strip()
bae = [p[len(PREFIX):] for p in git("ls-tree", "-r", "--name-only", "HEAD", PREFIX + "bae_q_20261024/core", PREFIX + "bae_q_20261024/source").decode().split("\n") if p]
extra = ["bae_q_20261024/inputs/independent/F1.mat",
         "round3_20261024/core/coherence_scores.m", "round3_20261024/core/apply_gain.m",
         "qscan_offset_20261025/scripts/make_dev.m", "qscan_offset_20261025/scripts/run_offset.m",
         "qscan_offset_20261025/scripts/summarize_offset.m", "qscan_offset_20261025/scripts/qscan_select.m",
         "qscan_offset_20261025/scripts/noise_replica.m", "qscan_offset_20261025/scripts/estimate_sigmas.m",
         "qscan_offset_20261025/analysis/offset_summary.csv"]
files = {}
for path in bae + extra:
    data = blob(path); out = DEST / path; out.parent.mkdir(parents=True, exist_ok=True); out.write_bytes(data)
    files[path] = dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                             commit=head, files=files), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files from {head[:8]}")
