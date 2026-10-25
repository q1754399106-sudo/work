"""Copy every input this round reads (including the previous round's sealed predictions and truth tables) into the ignored cache/ and verify it byte for byte.

Sources are Git objects of the current HEAD (regular blobs, or local Git LFS objects checked
against their pointer oid), so the run does not depend on the shared sparse checkout, which
other sessions rewrite. Writes cache/SNAPSHOT.json with commit, path, size and SHA-256.
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
bae = [p for p in git("ls-tree", "-r", "--name-only", "HEAD", "bae_q_20261024/core", "bae_q_20261024/source",
                      "bae_q_20261024/inputs/independent", "bae_q_20261024/truth/independent").decode().split()]
extra = ["qscan_criterion_20261024/inputs/confirm/F3.mat", "qscan_criterion_20261024/inputs/confirm/F4.mat",
         "qscan_criterion_20261024/truth/confirm/F3_truth.mat", "qscan_criterion_20261024/truth/confirm/F4_truth.mat",
         "round3_20261024/dev/F3.mat", "round3_20261024/dev/F4.mat",
         "round3_20261024/dev/F3_truth.mat", "round3_20261024/dev/F4_truth.mat",
         "qscan_bias_variance_20261024/analysis/decomposition.csv", "qscan_bias_variance_20261024/analysis/replicates.csv",
         "qscan_bias_variance_20261024/analysis/snr_probe.csv", "qscan_bias_variance_20261024/analysis/qscan_vs_base.csv",
         "qscan_bias_diagnostic_20261024/predictions/manifest.csv", "qscan_bias_diagnostic_20261024/analysis/expected_truth.csv",
         "qscan_bias_diagnostic_20261024/analysis/scored_gathers.csv", "bae_q_20261024/make_independent.m",
         *[f"qscan_bias_diagnostic_20261024/predictions/diagnostics_w{w}.csv" for w in range(1, 5)]]
files = {}
for path in bae + extra:
    data = blob(path)
    out = DEST / path
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    files[path] = dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                             commit=head, files=files), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files from {head[:8]}")
