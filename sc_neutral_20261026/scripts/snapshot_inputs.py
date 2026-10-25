"""Copy every input of PLAN_ZH.md into the ignored cache/ and verify it.

Git inputs come from the registered commits (LFS objects checked against their pointer oid) and
must equal the HEAD blob. The fd2d sweep outputs are local-only (fd2d_20261025/cache/sweep, never
uploaded); their SHA-256 is recorded so that the run is tied to these exact files.
Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DEST = ROOT / "cache" / "inputs"
SOURCES = {
    "c09c39a": ["fd2d_20261025/models", "fd2d_20261025/truth",
                "fd2d_20261025/scripts/noise_replica.m", "fd2d_20261025/scripts/estimate_sigmas.m"],
    "b59a6c3": ["bae_q_20261024/source"],
    "ff95a69": ["novelty_evidence_20261025/scripts/benchmark.py", "novelty_evidence_20261025/scripts/common.py",
                "novelty_evidence_20261025/analysis/selection.json", "novelty_evidence_20261025/analysis/environment.json"],
}
SWEEP = REPO / "fd2d_20261025" / "cache" / "sweep"
MODELS = ["H1", "H2", "H3", "H4", "H5", "H6"]
SCENES = ["Q100_14dB", "Q100_08dB", "Q070_14dB"]


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=True).stdout


def blob(rev, path):
    data = git("show", f"{rev}:{path}")
    if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
        f = dict(line.split(" ", 1) for line in data.decode().strip().splitlines())
        oid, size = f["oid"].split(":", 1)[1], int(f["size"])
        data = (REPO / ".git" / "lfs" / "objects" / oid[:2] / oid[2:4] / oid).read_bytes()
        assert len(data) == size and hashlib.sha256(data).hexdigest() == oid, path
    return data


files = {}
for rev, paths in SOURCES.items():
    for path in git("ls-tree", "-r", "--name-only", rev, *paths).decode().split():
        data = blob(rev, path)
        assert data == blob("HEAD", path), f"{path} changed after {rev}"
        out = DEST / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[path] = dict(commit=rev, size=len(data), sha256=hashlib.sha256(data).hexdigest())
sweep = {}
for m in MODELS:
    for s in SCENES:
        for r in range(7):
            p = SWEEP / f"{m}_{s}" / f"r{r:03d}.mat"
            sweep[p.relative_to(REPO).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), files=files, fd2d_sweep_sha256=sweep), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} Git files, {len(sweep)} sweep files")
