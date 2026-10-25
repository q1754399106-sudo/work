"""Copy the Git inputs of PLAN_ZH.md into the ignored cache/ and record the local-only inputs.

Git files are read at their registered commit (old path, before the c625fd6 reorganisation) and
must equal the HEAD blob at the new path under 地震衰减补偿研究/. Local-only inputs (fd2d sweep
cache, sc_neutral rebuilt gathers and SC outputs) are recorded by SHA-256.
Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent                      # 地震衰减补偿研究/
REPO = RESEARCH.parent
DEST = ROOT / "cache" / "inputs"
SOURCES = {
    "c09c39a": ["fd2d_20261025/models", "fd2d_20261025/truth"],
    "ff95a69": ["novelty_evidence_20261025/scripts/benchmark.py", "novelty_evidence_20261025/scripts/common.py",
                "novelty_evidence_20261025/analysis/selection.json"],
}
MODELS = [f"H{m}" for m in range(1, 7)]
SCENES = {"Q100_14dB": 1, "Q100_08dB": 2, "Q070_14dB": 1}


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


sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
files = {}
for rev, paths in SOURCES.items():
    for path in git("ls-tree", "-r", "--name-only", rev, *paths).decode().split():
        data = blob(rev, path)
        assert data == blob("HEAD", f"地震衰减补偿研究/{path}"), f"{path} changed after {rev}"
        out = DEST / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[path] = dict(commit=rev, size=len(data), sha256=hashlib.sha256(data).hexdigest())
local = {}
for m in MODELS:
    for s, level in SCENES.items():
        for r in range(7):
            for p in [RESEARCH / "fd2d_20261025/cache/sweep" / f"{m}_{s}" / f"r{r:03d}.mat",
                      RESEARCH / "sc_neutral_20261026/cache/noisy" / f"{m}_L{level}_r{r}.mat",
                      RESEARCH / "sc_neutral_20261026/cache/solved" / f"{m}_{s}_r{r}_SC.npz",
                      RESEARCH / "sc_neutral_20261026/cache/solved" / f"{m}_{s}_r{r}_SC.json"]:
                local[p.relative_to(RESEARCH).as_posix()] = sha(p)
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), files=files, local_sha256=local), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} Git files, {len(local)} local files")
