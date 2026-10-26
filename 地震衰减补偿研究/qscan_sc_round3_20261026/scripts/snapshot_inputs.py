"""Copy the Git inputs of PLAN_ZH.md into the ignored cache/ and record the local-only inputs.

The novelty SC code is read at its registered commit ff95a69 (old path, before the c625fd6
reorganisation) and must equal the HEAD blob at the new path under 地震衰减补偿研究/. round3 data,
archived results and logs are read at HEAD. The local novelty SC outputs at Q+0 (used only for the
consistency check) are recorded by SHA-256. Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT.parent                      # 地震衰减补偿研究/
REPO = RESEARCH.parent
DEST = ROOT / "cache" / "inputs"
SETS = ["F3", "F4", "V4", "V5"]
NOISES = ["Moderate", "Strong"]
ERRORS = [-50, -30, -10, 0, 10, 30, 50]
SOURCES = {
    "ff95a69": ["novelty_evidence_20261025/scripts/benchmark.py", "novelty_evidence_20261025/scripts/common.py",
                "novelty_evidence_20261025/analysis/selection.json"],
    "HEAD": [f"round3_20261024/dev/{s}{x}.mat" for s in SETS for x in ("", "_truth")]
            + [f"round3_20261024/results/{s}_Q{e:+d}_{n}.mat" for s in SETS for e in ERRORS for n in NOISES]
            + [f"round3_20261024/logs/{f}" for f in ("combo_F34.log", "combo_V4F.log", "combo_V5.log")],
}


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
    for path in paths:
        data = blob(rev, path if rev != "HEAD" else f"地震衰减补偿研究/{path}")
        if rev != "HEAD":
            assert data == blob("HEAD", f"地震衰减补偿研究/{path}"), f"{path} changed after {rev}"
        out = DEST / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[path] = dict(commit=rev, size=len(data), sha256=hashlib.sha256(data).hexdigest())
local = {}
for s in SETS:
    for n in NOISES:
        for ext in ("npz", "json"):
            p = RESEARCH / "novelty_evidence_20261025/cache/solved" / f"{s}_{n}_0.125_10.{ext}"
            local[p.relative_to(RESEARCH).as_posix()] = sha(p)
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), files=files, local_sha256=local), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} Git files, {len(local)} local files")
