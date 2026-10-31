"""PRE_RUN_HASHES.txt: SHA-256 of the plan, CHOICE.json, the experiment-17 scripts, and every sealed data file."""
import hashlib
from lib17 import ROOT, DATA3, TEST

files = [ROOT / "PLAN_ZH.md", ROOT / "CHOICE.json", *sorted((ROOT / "scripts").glob("*.py")), DATA3 / "analysis" / "datasets.json",
         *sorted((DATA3 / "scripts").glob("*.py")), *[DATA3 / "cache" / "data" / f"{d}.npz" for d in TEST]]
lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT.parent).as_posix()}" for p in files]
(ROOT / "PRE_RUN_HASHES.txt").write_bytes(("\n".join(lines) + "\n").encode("utf-8")); print(len(lines), "files")
