"""analysis/summary.csv and a markdown table of the data-level quantities (no method is run)."""
import csv, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; ANA = ROOT / "analysis"
meta = json.loads((ANA / "datasets.json").read_text(encoding="utf-8"))
real = json.loads((ANA / "realism.json").read_text(encoding="utf-8"))
rows = []
for did, m in sorted(meta.items(), key=lambda kv: (kv[1]["structure"], kv[1]["c0"], kv[1]["r"], kv[1]["variant"] != "Z")):
    lo = hi = ""
    if m["variant"] == "Z":
        lo, hi = (f"{v:.2f}" for v in np.percentile(m["c_eff_dif_trace"], [5, 95]))
    rows.append(dict(id=did, c0=m["c0"], c_eff_cum=f"{m['c_eff_cum']:.3f}", c_eff_dif=f"{m['c_eff_dif']:.3f}",
                     c_dif_trace_p5=lo, c_dif_trace_p95=hi, level=f"{m['level']:.2f}", f_true=f"{m['f_true']:.0f}",
                     decay_slope=f"{m['decay_slope']:.4f}", decay_slope_real=f"{real[m['structure']]['decay_slope_real']:.4f}",
                     snr_cal_err=f"{m.get('snr_cal_err', float('nan')):.3f}",
                     gas_shadow=("" if m["gas_shadow_30_60Hz"] is None else f"{m['gas_shadow_30_60Hz']:.2f}")))
with open(ANA / "summary.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("| 编号 | c₀ | c_eff 累积 | c_eff 增量 | 各道增量 5–95% | 噪声倍数 | f_true (Hz) | 衰减斜率（实测） | 气藏下方振幅比 |")
print("|---|---|---|---|---|---|---|---|---|")
for r in rows:
    span = f"{r['c_dif_trace_p5']}–{r['c_dif_trace_p95']}" if r["c_dif_trace_p5"] else "—"
    print(f"| {r['id']} | {r['c0']} | {r['c_eff_cum']} | {r['c_eff_dif']} | {span} | {r['level']} | {r['f_true']} | "
          f"{r['decay_slope']}（{r['decay_slope_real']}） | {r['gas_shadow'] or '—'} |")
