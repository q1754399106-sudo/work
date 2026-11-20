"""Line selection of the sixth batch (DATA_ZH.md; uqw_snr_20261114/PLAN_CONFIRM_S_ZH.md A): select_lines5.py of the fifth batch
with two changes only: (1) the fifth batch's candidate list lines5.json counts as a field-line list (it lists every candidate,
not only the chosen lines), and its chosen structure and noise lines count as used; (2) the fifth batch's structures join the
used structure lines. Same choice rule: structures = first, middle and last eligible candidates (>= 12 inline numbers from every
used structure line); noise lines = 6 of the remaining candidates >= 4 inline numbers from the chosen structures, evenly spaced.
Mentions in this batch's folder and in the confirmation plan are ignored."""
import subprocess, json, sys

P5 = "地震衰减补偿研究/semisynth_fault_q5_20261109/lines5.json"
FIELD_LISTS = {"地震衰减补偿研究/field_baselines_volve_20261029/analysis/lines.json",
               "地震衰减补偿研究/sca_ceiling_20261030/analysis/field_bands.json",
               "地震衰减补偿研究/field_vsp_volve_f11_20261026/probe/cube_coverage.json",
               "地震衰减补偿研究/qscan_anchor_theory_20261030/analysis/wlen.json", P5}
SELF = (":!地震衰减补偿研究/semisynth_fault_q6_20261114/*", ":!地震衰减补偿研究/uqw_snr_20261114/PLAN_CONFIRM_S_ZH.md")
USED_STRUCT = [10141, 10227, 10243, 10253, 10271, 10291, 10305, 10325, 10345, 10351, 10363, 10261, 10279, 10333, 10007, 10105, 10213]
L5 = json.loads(subprocess.run(["git", "-c", "core.quotepath=off", "show", "HEAD:" + P5], capture_output=True, text=True, encoding="utf-8").stdout)
USED5 = set(L5["structure"] + L5["noise"])


def mentions(L):
    out = subprocess.run(["git", "-c", "core.quotepath=off", "grep", "-l", "-w", str(L), "HEAD", "--", "*.py", "*.md", "*.json",
                          ":!*/results*/*", ":!*/cache/*", *SELF], capture_output=True, text=True, encoding="utf-8").stdout
    return {x.split(":", 1)[1] for x in out.splitlines() if x}


cand = [L for L in range(9999, 10220, 2) if L not in USED5 and mentions(L) <= FIELD_LISTS]
struct_ok = [L for L in cand if min(abs(L - u) for u in USED_STRUCT) >= 12]
struct = [struct_ok[0], struct_ok[len(struct_ok) // 2], struct_ok[-1]]
rest = [L for L in cand if min(abs(L - s) for s in struct) >= 4]
noise = [rest[round(i * (len(rest) - 1) / 5)] for i in range(6)]
res = dict(candidates=cand, structure_eligible=struct_ok, structure=struct, noise=noise)
print(json.dumps(res))
if "--write" in sys.argv:
    json.dump(res, open(sys.argv[sys.argv.index("--write") + 1], "w", encoding="utf-8"), indent=1)
