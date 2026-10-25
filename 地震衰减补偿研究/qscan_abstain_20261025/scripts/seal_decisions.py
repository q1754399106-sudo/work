"""Seal all truth-free decisions and saved outputs before score.m reads any target."""
from pathlib import Path
import datetime, hashlib, json

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "decisions"
assert not (D / "SEALED_SHA256.json").exists() and not (ROOT / "analysis" / "gather_scores.csv").exists()
csvs = sorted(D.glob("K*.csv"))
assert len(csvs) == 32, f"{len(csvs)} of 32 model x scene decision files"
outs = sorted((ROOT / "cache" / "out").glob("K*/r*.mat"))
assert len(outs) == 192, f"{len(outs)} of 192 outputs"
paths = [ROOT / "PLAN_ZH.md", *csvs, *sorted((ROOT / "scripts").glob("*")), *outs, *sorted((ROOT / "models").glob("*.mat"))]
rec = dict(sealed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), truth_read_before_seal=False,
           hashes={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(D / "SEALED_SHA256.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
print(f"sealed {len(paths)} files")
