"""Restore every input of experiment 12 into the ignored cache/inputs and verify it.

Same sources as experiment 11 (snapshot_inputs.py of srsc_volve_20261028) plus the experiment-11 code (sc11.py,
run11.py) and its development records. Git files come from their registered commits (LFS resolved through
`git lfs smudge`, checked against the pointer oid) and must equal the HEAD blob. Local files (MIG_RAW cube, the
experiment-10 held-out outputs H_*.npz and the experiment-11 development outputs cache/dev/H_*.npz, all seen data)
are recorded by SHA-256. Writes cache/SNAPSHOT.json (read by semi.read_lines for the cube path).
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
PREFIX = "地震衰减补偿研究/"
DEST = ROOT / "cache" / "inputs"
E6 = "field_vsp_volve_f11_20261026/"
E8 = "semisynthetic_volve_20261028/"
E9 = "multiwell_tie_volve_20261028/"
E10 = "specslope_q_volve_20261028/"
E11 = "srsc_volve_20261028/"
WELLS = ["F-1", "F-1_A", "F-1_B", "F-10", "F-11_A", "F-11_T2", "F-12", "F-14", "F-15", "F-15_A", "F-15_B", "F-4", "F-5"]
DEV = [f"H_c{c}_s{s}" for c in (0.55, 0.7, 0.9, 1.2, 1.4) for s in (9101, 9102, 9103, 9104)]
SOURCES = {  # commit -> (paths at that commit, prefixed?)
    "ff95a69": (["novelty_evidence_20261025/scripts/benchmark.py", "novelty_evidence_20261025/scripts/common.py",
                 "novelty_evidence_20261025/analysis/selection.json"], False),
    "0b04ea9": ([E6 + "inputs/posthoc_curtain.mat", E6 + "analysis/posthoc_references.mat"], True),
    "b1d76f8": (["field_application_20261028/scripts/lib.py", "field_application_20261028/analysis/inline.json"],
                True),
    "c395b78": ([E8 + "scripts/semi.py", E8 + "analysis/build.json"], True),
    "13e2fa7": ([E9 + "inventory/wells.csv", E9 + "raw/picks/Well_picks_Volve_v1.dat"]
                + [E9 + f"inventory/logs/15_9-{w}.csv" for w in WELLS], True),
    "346b562": ([E10 + "scripts/sq.py", E10 + "scripts/gen.py", E10 + "scripts/run.py", E10 + "results/real/lines.json"]
                + [E10 + f"results/synth/{d}.json" for d in DEV], True),
    "fbad95e": ([E9 + "scripts/wq.py", E9 + "scripts/synth.py", E9 + "scripts/tie.py", E9 + "results/tie.json"], True),
    "f7ed1d2": ([E11 + "scripts/sc11.py", E11 + "scripts/run11.py"], True),
    "8061a9e": ([E11 + f"analysis/dev/{d}.json" for d in DEV], True),
}
LOCAL = {
    "MIG_RAW.segy": (ROOT.parent / E6 / "cache" / "seismic" / "MIG_RAW.segy",
                     "22dbd614ee34a8d2becf3ea66d55db6887b8ed7804212452077aef9340c9ef92"),
}


def git(*args, stdin=None):
    return subprocess.run(["git", "-C", str(REPO), *args], input=stdin, capture_output=True, check=True).stdout


def blob(rev, path):
    data = git("show", f"{rev}:{path}")
    if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
        oid = dict(line.split(" ", 1) for line in data.decode().strip().splitlines())["oid"].split(":", 1)[1]
        data = git("lfs", "smudge", "--", path, stdin=data)
        assert hashlib.sha256(data).hexdigest() == oid, path
    return data


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


files = {}
for rev, (paths, prefixed) in SOURCES.items():
    for path in paths:
        data = blob(rev, PREFIX + path if prefixed else path)
        assert data == blob("HEAD", PREFIX + path), f"{path} changed after {rev}"
        out = DEST / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[path] = dict(commit=rev, size=len(data), sha256=hashlib.sha256(data).hexdigest())
local = {}
for k, (p, h) in LOCAL.items():
    assert sha(p) == h, f"{k} differs"
    local[k] = dict(path=str(p), sha256=h)
for d in DEV:
    for tag, p in (("E10", ROOT.parent / E10 / "cache" / "out" / f"{d}.npz"),
                   ("E11", ROOT.parent / E11 / "cache" / "dev" / f"{d}.npz")):
        local[f"{tag}/{d}.npz"] = dict(path=str(p), sha256=sha(p))
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), cube=str(LOCAL["MIG_RAW.segy"][0]), local=local, files=files),
    indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files, {len(local)} local")
