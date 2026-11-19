"""Line selection of the sixth semi-synthetic batch and of the field no-well lines (PLAN_CONFIRM_SHAPING_ZH.md section 3).
Reads only the text of the repository at HEAD (one git grep pass); no seismic data, no compensation, no well tie.

Rules (fixed in the plan before running):
  candidates  odd inlines 9985-10369 (the Volve MIG_RAW cube) whose every mention in the repository (.py/.md/.json at HEAD,
              results and cache folders excluded) is in an inline-list file of LIST_FILES; mentions in this folder ignored.
  used        structure and noise lines of batches 1-5 (USED_STRUCT, USED_NOISE) are never candidates.
  structure   6 lines among the candidates inside the experiment-13 range 9999-10219 that are >= 12 inline numbers from
              every used structure line, evenly spread over that eligible list (asserted >= 12 from each other).
  noise       12 lines among the remaining candidates inside 9999-10219, >= 4 inline numbers from the chosen structures,
              evenly spaced over that list; noise lines 2k / 2k+1 are realizations r = 0 / 1 of structure k (as fq.NOISE).
  field_nowell candidates outside 9999-10219 (never compensated in a field test, never a structure or noise line).
usage: select_lines6.py [--write PATH]
"""
import subprocess, json, sys, re
from collections import defaultdict

LIST_FILES = {"地震衰减补偿研究/field_baselines_volve_20261029/analysis/lines.json",
              "地震衰减补偿研究/sca_ceiling_20261030/analysis/field_bands.json",
              "地震衰减补偿研究/field_vsp_volve_f11_20261026/probe/cube_coverage.json",
              "地震衰减补偿研究/qscan_anchor_theory_20261030/analysis/wlen.json",
              "地震衰减补偿研究/semisynth_fault_q5_20261109/lines5.json",
              "地震衰减补偿研究/field_qperturb_20261113/analysis/ablation.json"}   # list of the field lines of the I2 ablation
SELF = "地震衰减补偿研究/uqw_shaping_confirm_20261119"
USED_STRUCT = [10141, 10227, 10243, 10253, 10271, 10291, 10305, 10325, 10345, 10351, 10363, 10261, 10279, 10333, 10007, 10105, 10213]
USED_NOISE = [9987, 10311, 9993, 10327, 10339, 10361, 9991, 10259, 9997, 10287, 10331, 10365, 10223, 10281, 10229, 10299, 10313,
              10337, 10249, 10265, 10285, 10293, 10343, 10357, 10019, 10053, 10083, 10127, 10169, 10205]


def all_mentions():
    """{inline: set of files} for every whole-word 99xx/10xxx number in the scanned files."""
    out = subprocess.run(["git", "-c", "core.quotepath=off", "grep", "-o", "-w", "-E", r"(99[89][0-9]|10[0-3][0-9][0-9])", "HEAD", "--",
                          "*.py", "*.md", "*.json", ":!*/results*/*", ":!*/cache/*", ":!" + SELF + "/*"],
                         capture_output=True, text=True, encoding="utf-8").stdout
    m = defaultdict(set)
    for x in out.splitlines():
        p = x.split(":")                            # HEAD:path:match
        if len(p) >= 3 and re.fullmatch(r"\d+", p[-1]):
            m[int(p[-1])].add(":".join(p[1:-1]))
    return m


M = all_mentions(); used = set(USED_STRUCT) | set(USED_NOISE)
cand = [L for L in range(9985, 10370, 2) if L not in used and M.get(L, set()) <= LIST_FILES]
inside = [L for L in cand if 9999 <= L <= 10219]
elig = [L for L in inside if min(abs(L - u) for u in USED_STRUCT) >= 12]
struct = [elig[round(i * (len(elig) - 1) / 5)] for i in range(6)]          # evenly spread over the eligible list
assert all(abs(a - b) >= 12 for i, a in enumerate(struct) for b in struct[i + 1:]), struct
rest = [L for L in inside if L not in struct and min(abs(L - s) for s in struct) >= 4]
noise = [rest[round(i * (len(rest) - 1) / 11)] for i in range(12)] if len(rest) >= 12 else []
res = dict(candidates=cand, inside=inside, structure_eligible=elig, structure=struct, noise=noise,
           field_nowell=[L for L in cand if not 9999 <= L <= 10219],
           pairs=[dict(structure=s, noise_r0=noise[2 * k] if noise else None, noise_r1=noise[2 * k + 1] if noise else None)
                  for k, s in enumerate(struct)])
print(json.dumps(res))
if "--write" in sys.argv:
    json.dump(res, open(sys.argv[sys.argv.index("--write") + 1], "w", encoding="utf-8"), indent=1)
