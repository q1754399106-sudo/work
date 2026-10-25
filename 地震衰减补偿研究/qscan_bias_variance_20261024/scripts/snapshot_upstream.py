"""Byte-copy the upstream solver code and F1/F2 data into the ignored cache/ and verify SHA-256.

The repository uses a sparse checkout that other jobs may narrow while MATLAB runs; the snapshot
keeps the run independent of that. bv_setup.m uses it only when SNAPSHOT.json is present.
"""
from pathlib import Path
import datetime, hashlib, json, shutil

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DEST = ROOT / "cache" / "bae_q_20261024"
recorded = json.loads((ROOT / "UPSTREAM_HASHES.json").read_text(encoding="utf-8"))["hashes"]
keep = ("bae_q_20261024/core/", "bae_q_20261024/source/", "bae_q_20261024/inputs/independent/",
        "bae_q_20261024/truth/independent/", "bae_q_20261024/results/ind/", "bae_q_20261024/results/ind_qscan/")
copied = {}
for name, digest in recorded.items():
    if not name.startswith(keep):
        continue
    src, dst = REPO / name, DEST / Path(name).relative_to("bae_q_20261024")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    assert hashlib.sha256(dst.read_bytes()).hexdigest() == digest, name
    copied[name] = digest
(DEST / "SNAPSHOT.json").write_text(json.dumps(dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                    files=copied), indent=2), encoding="utf-8")
print(f"snapshot verified: {len(copied)} files")
