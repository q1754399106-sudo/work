"""Well ties of the field Q-perturbation test (PLAN_ZH.md section 4): tie13.py of experiment 13 with the sections of
cache/out/IL_<L>.npz. Writes results/tie.json; prints no scores.
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run as R  # noqa: E402
from run import L, np, semi, lib, TW  # noqa: E402
import run13  # noqa: E402
from lib13 import tie  # noqa: E402


def main():
    logs = run13.well_logs(); G, mid = tie.groups(logs); h, x, y, il, xl = L.cube(); cm = semi.curtain_cm(); fp = float(cm["fpeak"]); sec = {}; rows = {}

    def section(ln):
        if ln not in sec:
            sec[ln] = dict(np.load(R.CACHE / f"IL_{ln}.npz"))
        return sec[ln]
    for n, g in logs.items():
        w, idx = tie.path(g, x, y, il, xl); rf = np.zeros(len(semi.TF)); tm = np.round(0.5 * (g["t"][1:] + g["t"][:-1]) - 0.0005, 6)
        i0 = int(round((tm[0] - semi.TF[0]) / L.wq.FINE)); rf[i0:i0 + len(g["rc"])] = g["rc"]
        ref = semi.ricker_ref(rf, fp); ws = semi.win(TW, g["t0"], g["t1"]); sc = {}
        for m in list(section(int(il[idx[0, 0]]))):
            tr = np.zeros(len(TW))
            for ii, r in zip(np.where(w)[0], idx):
                tr[ii] = np.mean([section(int(il[j]))[m][ii, (xl[j] - L.XL0) // L.XSTEP] for j in r])
            sc[m] = dict(zip(("T2", "shift_ms"), lib.tie(tr, ref, TW, ws)))
        rows[n] = dict(span=[g["t0"], g["t1"]], scores=sc); print(n, "done", flush=True)
    (R.ROOT / "results" / "tie.json").write_text(json.dumps(dict(groups=G, wells=rows), indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
