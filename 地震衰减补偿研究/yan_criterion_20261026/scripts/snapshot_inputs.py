"""Copy the Git inputs of PLAN_ZH.md (from HEAD) into the ignored cache/ and record the local-only
sweep caches of qscan_offset_20261025 and fd2d_20261025 by SHA-256. Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent
REPO = RESEARCH.parent
DEST = ROOT / "cache" / "inputs"
GIT = ["bae_q_20261024/source",
       "qscan_offset_20261025/models", "qscan_offset_20261025/truth",
       "qscan_offset_20261025/scripts/noise_replica.m", "qscan_offset_20261025/scripts/estimate_sigmas.m",
       "qscan_offset_20261025/scripts/qscan_select.m", "qscan_offset_20261025/analysis/offset_summary.csv",
       "qscan_offset_20261025/analysis/offset_curves.csv", "fd2d_20261025/analysis/fd_curves.csv",
       "fd2d_20261025/models", "fd2d_20261025/truth",
       "fd2d_20261025/scripts/noise_replica.m", "fd2d_20261025/scripts/estimate_sigmas.m",
       "fd2d_20261025/analysis/fd_summary.csv"]
SWEEPS = {"G": ("qscan_offset_20261025", 8), "H": ("fd2d_20261025", 6)}
SCENES = ["Q100_14dB", "Q100_08dB", "Q070_14dB"]


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=True).stdout


def blob(path):
    data = git("show", f"HEAD:{path}")
    if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
        f = dict(line.split(" ", 1) for line in data.decode().strip().splitlines())
        oid, size = f["oid"].split(":", 1)[1], int(f["size"])
        data = (REPO / ".git" / "lfs" / "objects" / oid[:2] / oid[2:4] / oid).read_bytes()
        assert len(data) == size and hashlib.sha256(data).hexdigest() == oid, path
    return data


head = git("rev-parse", "HEAD").decode().strip()
files = {}
for rel in GIT:
    for path in git("ls-tree", "-r", "--name-only", "HEAD", f"地震衰减补偿研究/{rel}").decode("utf-8").split("\n"):
        if not path:
            continue
        data = blob(path)
        short = path.split("/", 1)[1]
        out = DEST / short
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[short] = dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
local = {}
for g, (d, n) in SWEEPS.items():
    for m in range(1, n + 1):
        for s in SCENES:
            for r in range(7):
                p = RESEARCH / d / "cache" / "sweep" / f"{g}{m}_{s}" / f"r{r:03d}.mat"
                local[p.relative_to(RESEARCH).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), head=head, files=files,
    local_sha256=local), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} Git files at {head[:7]}, {len(local)} local sweep files")
