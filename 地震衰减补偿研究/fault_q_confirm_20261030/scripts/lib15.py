"""Shared code of experiment 15 (PLAN_ZH.md): the experiment-12/13 methods, with their frozen parameters, on the
new semi-synthetic data with faults and complex Q (semisynth_fault_q_20261030).

Imported unchanged from the experiment-13 input snapshot: lib12 (iq, tvw), sc11 (SC-a, QScan, M1) and everything
they import (experiments 7-11). New code: data loading, the method set, local window scores, theory quantities.
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT.parent
E13 = BASE / "field_baselines_volve_20261029"; IN = E13 / "cache" / "inputs"
DATA = BASE / "semisynth_fault_q_20261030"
sys.path.insert(0, str(BASE / "novelty_evidence_20261025" / "cache" / "deps"))   # scipy
for p in ("novelty_evidence_20261025/scripts", "field_application_20261028/scripts", "semisynthetic_volve_20261028/scripts",
          "specslope_q_volve_20261028/scripts", "multiwell_tie_volve_20261028/scripts", "srsc_volve_20261028/scripts",
          "baselines_volve_20261029/scripts"):
    sys.path.insert(0, str(IN / p))
import lib12, sc11, wq  # noqa: E402
from sc11 import np, lib, semi, sq, TW, DT  # noqa: E402

for mod in (lib12, sc11, sq, semi, wq):
    mod.ROOT, mod.IN = E13, IN
wq.E9 = IN / "multiwell_tie_volve_20261028"
sc11.E10 = lib12.E10 = IN / "specslope_q_volve_20261028"

SEL = json.loads((IN / "baselines_volve_20261029" / "analysis" / "dev_selection.json").read_text(encoding="utf-8"))
A_TUNED = SEL["B2"]["a"]                                   # 2^-13 (experiment 12)
G_IQ = SEL["B1"]["g_db"]                                   # 70 dB
META = json.loads((DATA / "analysis" / "datasets.json").read_text(encoding="utf-8"))
C0 = (0.75, 1.05, 1.4)
TEST = [f"{v}{k}_c{c}_r{r}" for v in "ZY" for k in "BC" for c in C0 for r in (0, 1)]
SMOKE = ["ZA_c1.05_r0"]


def load(did):
    return {k: v.astype(float) for k, v in np.load(DATA / "cache" / "data" / f"{did}.npz").items()}, META[did]


def methods(D, cm, Q, c_sr, c_qs, c_or, p):
    """PLAN_ZH.md section 3. c_sr: SR factor; c_qs: QScan factor; c_or: c_eff_cum of the truth (oracle, descriptive)."""
    sc = lambda c: lib.solve(cm, Q * c, p["wv"], D, p["sigma"], p["K"], p["live"], A_TUNED, lib.SEL_SC[1])[0]
    out = dict(Input=D, IQ30=lib12.iq(D, cm, Q, c_sr, c_sr, 30), IQ70=lib12.iq(D, cm, Q, c_sr, c_sr, G_IQ))
    out["SCa_SR"] = sc(c_sr)
    out["VDB"] = lib12.tvw(lib12.iq(D, cm, Q, None, c_sr, 0), p["wv"]["w_out"], SEL["B3"]["eps"], win=SEL["B3"]["win"])
    out["QSCa"] = sc11.sca(D, cm, Q, c_qs, p)
    out["QSCa_tuned"] = sc(c_qs)
    out["SCa_oracle"] = sc(c_or)
    out["IQ70_oracle"] = lib12.iq(D, cm, Q, c_or, c_or, G_IQ)
    return out


def score(o, ideal):
    s = sc11.m1(o, ideal); s["S"] = 0.5 * (s["M1"] + s["M1_lag"]); return s


def masks(arr, gas):
    """Local windows in 2.0-3.0 s (descriptive): within 3 traces of a fault plane / at least 15 traces away;
    0.05-0.40 s below the base of the gas lens on lens traces / same rows on traces outside the lens.
    gas is the Z twin's lens weight (Y has the same geometry but no lens)."""
    w = semi.win(TW, 2.0, 3.0)[:, None]; fd = arr["fault_dist"]
    m = dict(fault=(fd <= 3) & w, far=(fd >= 15) & w)
    cols = gas.max(axis=0) > 0.5; base = int(np.max(np.where(gas[:, cols].max(axis=1) > 0.5)[0]))
    rows = np.zeros(len(TW), bool); rows[base + 12:min(base + 101, len(TW))] = True
    m["below_gas"] = rows[:, None] & cols[None, :]; m["beside_gas"] = rows[:, None] & (gas.max(axis=0) < 0.01)[None, :]
    return m


def local(o, ideal, mk):
    return {k: float(np.corrcoef(o[v], ideal[v])[0, 1]) for k, v in mk.items()}


def theory(arr, cm):
    """Theory-weighted effective 1/c seen by anchored QScan (uniform phi2 on 1.2-3.0 s):
    u_QS = sum a (T - T(1.4)) / sum a (a - a(1.4)), T = trace-mean true t*; equals 1/c for T = a/c.
    kappa = a(1.4) sum a / sum a^2 (0.675 for the model Q)."""
    a = lib.absorption(cm, TW); i = int(np.argmin(np.abs(TW - 1.4))); T = arr["tstar"].mean(axis=1)
    return dict(u_qs=float(a @ (T - T[i]) / (a @ (a - a[i]))), kappa=float(a[i] * a.sum() / (a @ a)))
