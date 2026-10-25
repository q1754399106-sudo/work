"""Seal the 48 truth-free confirmation predictions before any confirmation truth is computed."""
from pathlib import Path
import csv, datetime, hashlib, json

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "predictions"
assert not (P / "SEALED_SHA256.json").exists(), "Already sealed."
ids = sorted(int(r["id"]) for f in sorted(P.glob("confirm_w*.csv")) for r in csv.DictReader(f.open(encoding="utf-8")))
assert ids == list(range(1, 49)), f"diagnosed {len(ids)} of 48 gathers"
paths = [ROOT / "PLAN_ZH.md", ROOT / "DEV_PLAN_ZH.md", ROOT / "analysis" / "dev_selection.json", *sorted(P.glob("confirm_w*.csv")),
         *sorted((ROOT / "inputs").glob("*.mat")), ROOT / "cache" / "SNAPSHOT.json",
         *sorted((ROOT / "scripts").glob("*.m")), *sorted((ROOT / "scripts").glob("*.py")),
         *sorted((ROOT / "cache" / "outputs").glob("c*.mat"))]
receipt = dict(sealed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), truth_read_before_seal=False,
               note="truth/*.mat exist since generation but are read only by run_truth after this seal.",
               hashes={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(P / "SEALED_SHA256.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
print(f"sealed {len(paths)} files, 48 gathers")
