"""Seal the truth-free predictions before any scoring: every manifest gather must be diagnosed once."""
from pathlib import Path
import csv, datetime, hashlib, json

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "predictions"
assert not (P / "SEALED_SHA256.json").exists(), "Already sealed."
ids = [int(r["id"]) for r in csv.DictReader((P / "manifest.csv").open(encoding="utf-8"))]
done = [int(r["id"]) for f in sorted(P.glob("diagnostics_w*.csv")) for r in csv.DictReader(f.open(encoding="utf-8"))]
assert sorted(done) == sorted(ids), f"diagnosed {len(done)} of {len(ids)} gathers"
paths = [ROOT / "PLAN_ZH.md", P / "manifest.csv", *sorted(P.glob("diagnostics_w*.csv")), ROOT / "cache" / "SNAPSHOT.json",
         *sorted((ROOT / "scripts").glob("*.m")), *sorted((ROOT / "scripts").glob("*.py")),
         *sorted((ROOT / "cache" / "outputs").glob("g*.mat"))]
receipt = dict(sealed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), truth_read_before_seal=False,
               note="Estimator outputs g*.mat stay in the local cache (rebuildable); their hashes are recorded here.",
               hashes={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(P / "SEALED_SHA256.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
print(f"sealed {len(paths)} files, {len(done)} gathers")
