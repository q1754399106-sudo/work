"""PRE_RUN_HASHES.txt: SHA-256 of the experiment-16 scripts, the plan, and every sealed data file."""
import hashlib
from lib16 import ROOT, DATA, TEST

files = [ROOT / "PLAN_ZH.md", *sorted((ROOT / "scripts").glob("*.py")), DATA / "analysis" / "datasets.json",
         *sorted((DATA / "scripts").glob("*.py")), *[DATA / "cache" / "data" / f"{d}.npz" for d in TEST]]
lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT.parent).as_posix()}" for p in files]
(ROOT / "PRE_RUN_HASHES.txt").write_bytes(("\n".join(lines) + "\n").encode("utf-8")); print(len(lines), "files")
