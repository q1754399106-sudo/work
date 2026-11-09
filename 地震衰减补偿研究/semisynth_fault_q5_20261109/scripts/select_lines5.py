"""Line selection of the fifth batch (DATA_ZH.md; qunknown_comp_20261109/PLAN_CONFIRM_W_ZH.md A).
Candidates: odd inlines 9999-10219 (experiment-13 field range) whose every mention in the repository (.py/.md/.json at
HEAD, results and cache folders excluded) is in one of the field-line lists FIELD_LISTS. Structure lines also keep
>= 12 inline numbers from every structure line used before. Choice: structures = first, middle and last eligible
structure candidates; noise lines = 6 of the remaining candidates at least 4 inline numbers from the chosen structures,
evenly spaced over the list. Mentions in the confirmation plan itself and in this batch's folder are ignored."""
import subprocess, json, sys

FIELD_LISTS = {"地震衰减补偿研究/field_baselines_volve_20261029/analysis/lines.json",
               "地震衰减补偿研究/sca_ceiling_20261030/analysis/field_bands.json",
               "地震衰减补偿研究/field_vsp_volve_f11_20261026/probe/cube_coverage.json",
               "地震衰减补偿研究/qscan_anchor_theory_20261030/analysis/wlen.json"}
SELF_PLAN = "地震衰减补偿研究/qunknown_comp_20261109/PLAN_CONFIRM_W_ZH.md"     # mentions in this plan and batch do not count
USED_STRUCT = [10141, 10227, 10243, 10253, 10271, 10291, 10305, 10325, 10345, 10351, 10363, 10261, 10279, 10333]


def mentions(L):
    out = subprocess.run(["git", "-c", "core.quotepath=off", "grep", "-l", "-w", str(L), "HEAD", "--", "*.py", "*.md", "*.json",
                          ":!*/results*/*", ":!*/cache/*", ":!" + SELF_PLAN, ":!地震衰减补偿研究/semisynth_fault_q5_20261109/*"], capture_output=True, text=True, encoding="utf-8").stdout
    return {x.split(":", 1)[1] for x in out.splitlines() if x}


cand = [L for L in range(9999, 10220, 2) if mentions(L) <= FIELD_LISTS]
struct_ok = [L for L in cand if min(abs(L - u) for u in USED_STRUCT) >= 12]
struct = [struct_ok[0], struct_ok[len(struct_ok) // 2], struct_ok[-1]]
rest = [L for L in cand if min(abs(L - s) for s in struct) >= 4]
noise = [rest[round(i * (len(rest) - 1) / 5)] for i in range(6)]
res = dict(candidates=cand, structure_eligible=struct_ok, structure=struct, noise=noise)
print(json.dumps(res))
if "--write" in sys.argv:
    json.dump(res, open(sys.argv[sys.argv.index("--write") + 1], "w", encoding="utf-8"), indent=1)
