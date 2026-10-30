"""PRE_RUN_HASHES.txt: SHA-256 of the experiment-15 scripts, the plan, and every sealed data file."""
import hashlib
from pathlib import Path
from lib15 import ROOT, DATA, TEST

files = [ROOT / "PLAN_ZH.md", *sorted((ROOT / "scripts").glob("*.py")),
         DATA / "analysis" / "datasets.json", *[DATA / "cache" / "data" / f"{d}.npz" for d in TEST]]
lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT.parent).as_posix()}" for p in files]
(ROOT / "PRE_RUN_HASHES.txt").write_text("\n".join(lines) + "\n", encoding="utf-8"); print(len(lines), "files")
