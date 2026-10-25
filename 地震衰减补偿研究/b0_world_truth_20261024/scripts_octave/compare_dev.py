"""Stage 0 on Octave: compare dev_oct_*.csv with the MATLAB development values (analysis/dev_w*.csv)."""
from pathlib import Path
import json
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
m = pd.concat([pd.read_csv(f) for f in sorted((ROOT / "analysis").glob("dev_w*.csv"))]).set_index("id").sort_index()
o = pd.concat([pd.read_csv(f) for f in sorted((ROOT / "octave_check").glob("dev_oct_*.csv"))]).set_index("id").sort_index()
assert o.index.is_unique and set(o.index) == set(m.index), (len(o), len(m))
tau = json.loads((ROOT / "analysis" / "dev_selection.json").read_text())["selected"]["threshold"]
cols = [c for c in o.columns if c in m.columns and c != "seconds"]
diff = {c: float(np.max(np.abs(o[c] - m.loc[o.index, c]))) for c in cols}
out = dict(gathers=len(o), max_abs_diff=diff, qscan_factor_identical=bool((o.qscan_factor.round(12) == m.loc[o.index].qscan_factor.round(12)).all()),
           B0M_calls_identical_at_tau=int(((o.B0M_rel_delta_mse > tau) == (m.loc[o.index].B0M_rel_delta_mse > tau)).sum()),
           B0_calls_identical_at_0=int(((o.B0_rel_delta_mse > 0) == (m.loc[o.index].B0_rel_delta_mse > 0)).sum()),
           B0C_calls_identical_at_0=int(((o.B0C_rel_delta_mse > 0) == (m.loc[o.index].B0C_rel_delta_mse > 0)).sum()),
           seconds_octave_mean=float(o.seconds.mean()), seconds_matlab_mean=float(m.loc[o.index].seconds.mean()))
(ROOT / "octave_check" / "DEV_COMPARISON.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
