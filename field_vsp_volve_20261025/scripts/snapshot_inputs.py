"""Copy the frozen code named in PLAN_ZH.md into the ignored cache/ and verify it byte for byte.

Each file is read from the Git object of its registered commit and must equal the HEAD blob,
so the run does not depend on the shared sparse checkout. Writes cache/SNAPSHOT.json.
"""
from pathlib import Path
import datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DEST = ROOT / "cache" / "inputs"
SOURCES = {  # registered commit -> paths (PLAN_ZH.md, "方法")
    "b59a6c3": ["bae_q_20261024/core", "bae_q_20261024/source"],
    "f38f3e3": ["round3_20261024/core/coherence_scores.m", "round3_20261024/core/apply_gain.m"],
    "16aabfe": ["combo_bias_variance_20261024/scripts/combo_outputs.m",
                "combo_bias_variance_20261024/scripts/qscan_select.m"],
}


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=True).stdout


files = {}
for commit, paths in SOURCES.items():
    for path in git("ls-tree", "-r", "--name-only", commit, *paths).decode().split():
        data = git("show", f"{commit}:{path}")
        assert data == git("show", f"HEAD:{path}"), f"{path} changed after {commit}"
        out = DEST / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(data)
        files[path] = dict(commit=commit, size=len(data), sha256=hashlib.sha256(data).hexdigest())
(ROOT / "cache" / "SNAPSHOT.json").write_text(json.dumps(dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    head=git("rev-parse", "HEAD").decode().strip(), files=files), indent=2), encoding="utf-8")
print(f"snapshot: {len(files)} files")
