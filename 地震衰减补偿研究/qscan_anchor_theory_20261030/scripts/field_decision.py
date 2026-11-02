"""Post hoc (FIELD_DECISION_ZH.md): does the two-anchor diagnostic change the Q choice at the two Volve VSP wells?
Inputs already archived: results9 (two-anchor scans on the experiment-13 lines), the experiment-10 SR factor
(lib13.c_sr), the robust VSP intervals (analysis/vsp_robust.json). Lines: 10087 (2 inlines from the F-15A VSP
curtain, inline 10085) and 10141 (inside the F-11 T2 RSVSP footprint, inlines 10139-10151). No new scans.
-> analysis/field_decision.json"""
import json, sys
sys.path.insert(0, "../../field_baselines_volve_20261029/scripts")
from env15 import np  # noqa: F401
import lib13 as L13

VSP = {"F-15A (Theil-Sen)": (10087, 1.42, 1.02, 3.87), "F-15A (OLS)": (10087, 1.34, 1.00, 2.46),
       "F-11 T2 RS only (deep from model)": (10141, 1.01, 0.81, 1.29), "F-11 T2 both segments": (10141, 0.84, 0.34, 3.69)}
out = {}
for name, (ln, cv, lo, hi) in VSP.items():
    v = {a: json.load(open(f"../results9/IL_{ln}__w1.4_a{a}.json"))["v_sel"] for a in (0.7, 1.0, 1.4)}
    k = (v[0.7] - v[1.4]) / (1 / 0.7 - 1 / 1.4); cq = 1 / v[1.0]; csr = float(L13.c_sr(ln)); cinv = (1 - k) / (v[1.0] - k)
    inside = lambda c: bool(lo <= c <= hi)
    out[name] = dict(inline=ln, vsp_c=cv, vsp_90=[lo, hi], kappa_hat=k, c_qscan=cq, c_qscan_anchor_range=[1 / v[0.7], 1 / v[1.4]],
                     c_sr=csr, c_one_step_inverse=cinv, qscan_inside=inside(cq), sr_inside=inside(csr), inverse_inside=inside(cinv),
                     decision="SR" if 1 - k < 0.5 else "QScan")
    print(f"{name:34s} IL{ln} kappa {k:.3f} QScan {cq:.3f}{'*' if inside(cq) else ' '} SR {csr:.3f}{'*' if inside(csr) else ' '} "
          f"inverse {cinv:.3f}{'*' if inside(cinv) else ' '} anchors {1/v[0.7]:.2f}-{1/v[1.4]:.2f}  VSP {cv} [{lo}, {hi}]")
json.dump(out, open("../analysis/field_decision.json", "w"), indent=1)
