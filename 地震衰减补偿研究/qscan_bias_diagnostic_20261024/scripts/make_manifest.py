"""Write predictions/manifest.csv: the 80 gathers diagnosed in this round (frozen with the plan)."""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
rows = []


def add(group, set_, model, condition, snr, q_err, source, seed):
    rows.append(dict(id=len(rows) + 1, group=group, set=set_, model=model, condition=condition, snr_db=snr,
                     q_err=q_err, source=source, noise_seed=seed, boot_seed=1_300_000_000 + 1000 * (len(rows) + 1)))


SNR = {"Moderate": 14, "Strong": 8}
for m, model in [(1, "F1"), (2, "F2")]:
    for k, cond in [(1, "Moderate"), (2, "Strong")]:
        for r in range(1, 4):   # seeds of qscan_bias_variance_20261024 replicates 1-3
            add("dev", "bae", model, cond, SNR[cond], 0, "replicate", 1_000_000_000 + 1_000_000 * m + 100_000 * k + r)
    for level, snr in [(1, 20), (2, 26)]:
        for r in range(1, 3):   # seeds of its high-SNR probe, replicates 1-2
            add("dev", "bae", model, f"SNR{snr}", snr, 0, "replicate", 1_000_000_000 + 1_000_000 * m + 100_000 * (2 + level) + r)
    for cond in ["Moderate", "Strong"]:
        for q in [-30, 30]:
            add("dev", "bae", model, cond, SNR[cond], q, "archived", 0)
for s, set_ in [(2, "qc"), (3, "r3")]:
    for m, model in [(3, "F3"), (4, "F4")]:
        for k, cond in [(1, "Moderate"), (2, "Strong")]:
            add("confirm", set_, model, cond, SNR[cond], 0, "archived", 0)
            for r in (1, 2):    # replicates 1-2 of the confirmation truth replicates (1-20)
                add("confirm", set_, model, cond, SNR[cond], 0, "replicate",
                    1_100_000_000 + 10_000_000 * s + 1_000_000 * m + 100_000 * k + r)
            for q in [-30, 30]:
                add("confirm", set_, model, cond, SNR[cond], q, "archived", 0)
with (ROOT / "predictions" / "manifest.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(len(rows), "gathers")
