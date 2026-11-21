"""POST HOC (after s2_source V0-V5): ideal noise-free model with the Q error in the t0-domain operator (exponent eps*a, not eps*a/m).
K1op: centroid form, sensitivity a * fbar(w_K1) (vs (a/m) fbar in s2_ideal.py). SHop: amplitude-domain Shuey fit, NMO with stretch mute;
aligned with the NMO mute ('mute') and with all offsets ('full', as gate A2's aligned scene). -> results/s2_ideal_operr.json"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT.parent / "offset_q_avo_20261118" / "scripts" / "legacy2b"))
import gate1 as g1  # noqa: E402
MODELS = {"M24": [42, 70, 55, 95], "M25": [75, 45, 30, 65], "M26": [35, 40, 85, 50], "M27": [95, 60, 45, 35]}
EV = np.round(g1.C.event_t / g1.DT).astype(int); F = np.arange(1, 1001) * 0.1
ric = lambda f, fp: f ** 2 * np.exp(-(f / fp) ** 2)  # noqa: E731
bp = lambda f, fp: np.clip((f - 5) / 5, 0, 1) * np.clip((fp + 8 - f) / 8, 0, 1)  # noqa: E731
out = []
for mn, Q in MODELS.items():
    a0, ax = g1.a0_ax(np.array(Q, float))
    for fp in (40.0, 60.0):
        R0, B2 = ric(F, fp), bp(F, fp) ** 2
        for j, k0 in enumerate(EV):
            a = a0[k0] + ax[k0]; s = g1.SIN2[k0]
            m = np.array([np.gradient(g1.GEO[i]["tx"], g1.DT)[k0] for i in range(g1.NX)])
            mute = np.array([g1.GEO[i]["tx"][k0] <= g1.MUTE * g1.T[k0] for i in range(g1.NX)]); full = np.ones(g1.NX, bool)
            co = g1.rpp(j, np.array([g1.GEO[i]["p"][k0] * g1.C.vp[j] for i in range(g1.NX)]))
            wk1 = [R0 * B2 * ric(F / m[i], fp) / m[i] * np.exp(-np.pi * F * a[i] * (1 / m[i] - 1)) for i in range(g1.NX)]
            cen = lambda w: (w * F).sum() / w.sum()  # noqa: E731
            def sl(y, keep):
                sc = s[keep] - s[keep].mean(); return float((y[keep] * sc).sum() / (sc ** 2).sum())
            s_al = np.array([a[i] * cen(R0 ** 2 * B2) for i in range(g1.NX)]); s_k1 = np.array([a[i] * cen(wk1[i]) for i in range(g1.NX)])
            def ga(eps, nmo, keep):
                y = []
                for i in range(g1.NX):
                    w = (wk1[i] if nmo else R0 ** 2 * B2) * np.exp(-np.pi * F * eps * a[i]); y.append(co[i] * w.sum() / (R0 ** 2 * B2).sum())
                X = np.c_[np.ones(keep.sum()), s[keep]]; A_, B_ = np.linalg.lstsq(X, np.array(y)[keep], rcond=None)[0]; return B_ / A_
            S = lambda nmo, keep: (ga(0.01, nmo, keep) - ga(-0.01, nmo, keep)) / 0.02  # noqa: E731
            r = dict(model=mn, fp=fp, event=j, phi_K1op=1 - sl(s_k1, mute) / sl(s_al, mute), phi_K1op_full=1 - sl(s_k1, mute) / sl(s_al, full),
                     SH_al_mute=S(0, mute), SH_al_full=S(0, full), SH_nmo=S(1, mute), fbar_ratio_far=float(cen(wk1[np.where(mute)[0][-1]]) / cen(R0 ** 2 * B2)), m_far=float(m[mute][-1]))
            r["phi_SH_mute"] = 1 - r["SH_nmo"] / r["SH_al_mute"]; r["phi_SH_full"] = 1 - r["SH_nmo"] / r["SH_al_full"]; out.append(r)
med = {k: float(np.median([r[k] for r in out])) for k in out[0] if k not in ("model", "fp", "event")}
(ROOT / "results" / "s2_ideal_operr.json").write_text(json.dumps(dict(rows=out, median=med), indent=1), encoding="utf-8"); print(json.dumps(med, indent=1))
