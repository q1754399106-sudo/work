"""Copy every input of this round from Git objects of HEAD into the ignored cache/ (byte-verified),
so the run does not depend on the working tree. Writes cache/SNAPSHOT.json."""
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
listed = git("ls-tree", "-r", "--name-only", "HEAD", "bae_q_20261024/core", "bae_q_20261024/source",
             "bae_q_20261024/inputs/independent", "bae_q_20261024/truth/independent",
             "qscan_offset_20261025/models").decode().split()
extra = ["qscan_offset_20261025/analysis/offset_curves.csv", "qscan_offset_20261025/analysis/offset_summary.csv",
         "combo_bias_variance_20261024/analysis/sweep_curves.csv", "combo_bias_variance_20261024/analysis/replicates.csv"]
files = {}
for path in listed + extra:
    data = blob(path)
    out = DEST / path
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    files[path] = dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                             commit=head, files=files), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files from {head[:8]}")
