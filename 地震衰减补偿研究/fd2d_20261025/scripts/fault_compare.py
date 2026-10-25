"""Near-fault vs far errors from analysis/fault_windows.csv (summarize_fault.m).

For each model, scene and quantity (D, B, A, noise-free B): the near(h)/far - 1 change, and where
near(h) falls among the 10-trace control blocks cNN (fraction of controls with a lower error, and
whether near(h) lies outside the control range). Also the error-optimal factor in near5 and far,
in grid steps. Writes analysis/fault_compare.csv and prints the summary used in FAULT_ZH.md.
"""
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STEP = math.log(math.exp((math.log(2.2) - math.log(0.6)) / 12))
QTY = {"mse_D": "D", "mse_B": "B", "mse_A": "A", "mse_B_noise_free": "B0"}

a = pd.read_csv(ROOT / "analysis" / "fault_windows.csv")
rows = []
for (model, scene), g in a.groupby(["model", "scene"], sort=False):
    w = g.set_index("window")
    ctrl = g[g.window.str.match(r"c\d\d")]
    for h in ("near3", "near5", "near10"):
        for col, name in QTY.items():
            v, far = w.at[h, col], w.at["far", col]
            rows.append({
                "model": model, "scene": scene, "window": h, "quantity": name,
                "near": v, "far": far, "rel_near_far": v / far - 1,
                "ctrl_min": ctrl[col].min(), "ctrl_max": ctrl[col].max(),
                "ctrl_frac_below": (ctrl[col] < v).mean(),
                "outside_ctrl": "above" if v > ctrl[col].max() else "below" if v < ctrl[col].min() else "",
                "n_ctrl": len(ctrl),
                "q_best_near": w.at[h, "q_best"], "q_best_far": w.at["far", "q_best"],
                "q_shift_steps": math.log(w.at[h, "q_best"] / w.at["far", "q_best"]) / STEP,
            })
C = pd.DataFrame(rows)
C.to_csv(ROOT / "analysis" / "fault_compare.csv", index=False, float_format="%.6g")

pd.set_option("display.width", 200)
n5 = C[C.window == "near5"]
print(n5.pivot_table(index=["scene", "model"], columns="quantity", values="rel_near_far").round(3)
      [["B0", "B", "D", "A"]].to_string())
print("\noutside control range (near5):")
print(n5[n5.outside_ctrl != ""][["scene", "model", "quantity", "rel_near_far", "outside_ctrl"]].to_string(index=False))
print("\ncontrol fraction below near5:")
print(n5.pivot_table(index=["scene", "model"], columns="quantity", values="ctrl_frac_below").round(2)
      [["B0", "B", "D", "A"]].to_string())
q = n5[n5.quantity == "B"][["scene", "model", "q_best_far", "q_best_near", "q_shift_steps"]]
print("\noptimum shift near5 vs far (grid steps):")
print(q.round(3).to_string(index=False))
