"""Well ties of experiment 13 (PLAN_ZH.md section 4), as experiment 9 (tie.py run) except that every trace is
taken from its compensated inline section (cache/out/IL_<L>_{A,B}.npz) instead of being compensated alone.

Per well: the log reflectivity and anchored span of tie.tie_log, the Ricker reference of tie.run, and at every 4 ms
sample the mean of the nearest trace and its four neighbours (inline +-2, xline +-2) along tie.path; score T2 =
lib.tie (Pearson over the span, best of +-20 ms). Writes results/tie.json; prints no scores.
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib13 as L  # noqa: E402
from lib13 import np, semi, lib, tie, TW  # noqa: E402
import run13  # noqa: E402

CACHE = L.ROOT / "cache" / "out"


def main():
    logs = run13.well_logs(); G, mid = tie.groups(logs); h, x, y, il, xl = L.cube()
    cm = semi.curtain_cm(); fp = float(cm["fpeak"]); sec = {}

    def section(ln):
        if ln not in sec:
            sec[ln] = {**dict(np.load(CACHE / f"IL_{ln}_A.npz")), **dict(np.load(CACHE / f"IL_{ln}_B.npz"))}
        return sec[ln]
    rows = {}
    for n, g in logs.items():
        w, idx = tie.path(g, x, y, il, xl)
        rf = np.zeros(len(semi.TF)); tm = np.round(0.5 * (g["t"][1:] + g["t"][:-1]) - 0.0005, 6)
        i0 = int(round((tm[0] - semi.TF[0]) / L.wq.FINE)); rf[i0:i0 + len(g["rc"])] = g["rc"]
        ref = semi.ricker_ref(rf, fp); ws = semi.win(TW, g["t0"], g["t1"])
        names = list(section(int(il[idx[0, 0]])))
        trs, sc = {}, {}
        for m in names:
            tr = np.zeros(len(TW))
            for ii, r in zip(np.where(w)[0], idx):
                tr[ii] = np.mean([section(int(il[j]))[m][ii, (xl[j] - L.XL0) // L.XSTEP] for j in r])
            trs[m] = tr; sc[m] = dict(zip(("T2", "shift_ms"), lib.tie(tr, ref, TW, ws)))
        jm = int(np.argmin((x - mid[n][0]) ** 2 + (y - mid[n][1]) ** 2))
        rows[n] = dict(span=[g["t0"], g["t1"]], inline_mid=int(il[jm]), n_traces=int(len(np.unique(idx))), scores=sc)
        np.savez_compressed(CACHE / f"tie_{n}.npz", ref=ref, **trs)
        print(n, "done", flush=True)
    (L.ROOT / "results").mkdir(exist_ok=True)
    (L.ROOT / "results" / "tie.json").write_text(json.dumps(dict(groups=G, seen=tie.SEEN, fpeak=fp, wells=rows), indent=1),
                                                encoding="utf-8")


if __name__ == "__main__":
    main()
