"""analysis/bands.json: shallow and target band edges (lib.wavelets and the robust variant) against f_true for the
18 Z datasets. Data-level only (no compensation is run)."""
import json
from libc import ROOT, IDS, META, load, wavelets, semi

cm = semi.curtain_cm(); rows = {}
for did in IDS:
    D = load(did)[0]["D"]; w0 = wavelets(D, cm, False); w1 = wavelets(D, cm, True)
    rows[did] = dict(f_true=META[did]["f_true"], noise_line=META[did]["noise_line"], level=META[did]["level"],
                     fs=w0["fmax_shallow"], fs_robust=w1["fmax_shallow"], fd=w0["fmax_target"], f_peak=w0["f_peak"])
    print(did, rows[did], flush=True)
(ROOT / "analysis").mkdir(exist_ok=True)
(ROOT / "analysis" / "bands.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
