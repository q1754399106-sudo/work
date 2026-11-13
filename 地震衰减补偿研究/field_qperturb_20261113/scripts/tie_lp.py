"""Post hoc, descriptive (not part of PLAN_ZH.md; verdicts R1-R5 are unchanged): the well ties of tie.py after a zero-phase
low-pass of every section (all methods and the input; pass <= 70 Hz, cosine taper to zero at 85 Hz). Question: how much of the
UQ-W / IQ70 deficit against VDB comes from the amplified energy above the signal band? Cut-off fixed before any tie was computed.
Writes results/tie_lp.json and prints the group statistics of analyze.py for comparison.
"""
import sys, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run as R  # noqa: E402
from run import L, np, semi, lib, TW, DT  # noqa: E402
import run13  # noqa: E402
from lib13 import tie  # noqa: E402

F0, F1 = 70.0, 85.0; M_ALL = R.M_ALL; LOW, HIGH = (0.5, 0.6, 0.7), (1.4, 1.7, 2.0)


def lowpass(x):
    n = x.shape[0]; nf = 2 * n; f = np.fft.rfftfreq(nf, DT); w = np.clip((F1 - f) / (F1 - F0), 0, 1); w = 0.5 - 0.5 * np.cos(np.pi * w)
    return np.fft.irfft(np.fft.rfft(x, nf, axis=0) * w[:, None], nf, axis=0)[:n]


def main():
    logs = run13.well_logs(); G, mid = tie.groups(logs); h, x, y, il, xl = L.cube(); cm = semi.curtain_cm(); fp = float(cm["fpeak"]); sec = {}; rows = {}

    def section(ln):
        if ln not in sec:
            sec[ln] = {k: lowpass(v.astype(float)) for k, v in np.load(R.CACHE / f"IL_{ln}.npz").items()}
        return sec[ln]
    for n, g in logs.items():
        w, idx = tie.path(g, x, y, il, xl); rf = np.zeros(len(semi.TF)); tm = np.round(0.5 * (g["t"][1:] + g["t"][:-1]) - 0.0005, 6)
        i0 = int(round((tm[0] - semi.TF[0]) / L.wq.FINE)); rf[i0:i0 + len(g["rc"])] = g["rc"]
        ref = semi.ricker_ref(rf, fp); ws = semi.win(TW, g["t0"], g["t1"]); sc = {}
        for m in list(section(int(il[idx[0, 0]]))):
            tr = np.zeros(len(TW))
            for ii, r in zip(np.where(w)[0], idx):
                tr[ii] = np.mean([section(int(il[j]))[m][ii, (xl[j] - L.XL0) // L.XSTEP] for j in r])
            sc[m] = lib.tie(tr, ref, TW, ws)[0]
        rows[n] = sc
    (R.ROOT / "results" / "tie_lp.json").write_text(json.dumps(dict(groups=G, wells=rows, lowpass=[F0, F1]), indent=1), encoding="utf-8")
    old = json.loads((R.ROOT / "results" / "tie.json").read_text(encoding="utf-8"))["wells"]; grp = lambda fn: np.array([np.mean([fn(n) for n in mem]) for mem in G.values()])

    def st(d):
        m = d.mean(); se = d.std(ddof=1) / np.sqrt(len(d)); return "%+.3f [%+.3f, %+.3f] %d/8" % (m, m - 1.895 * se, m + 1.895 * se, (d > 0).sum())
    for tag, t2 in (("original", lambda n, a, m: old[n]["scores"][f"{a}_{m}"]["T2"]), ("low-pass", lambda n, a, m: rows[n][f"{a}_{m}"])):
        print(tag, " input %.3f" % grp(lambda n: (old[n]["scores"]["Input"]["T2"] if tag == "original" else rows[n]["Input"])).mean())
        for a in R.METHODS:
            print("  %-5s T2 by m:" % a, " ".join("%.3f" % grp(lambda n: t2(n, a, m)).mean() for m in M_ALL))
        for name, ms, b in (("low  UQW-VDB ", LOW, "VDB"), ("low  UQW-IQ70", LOW, "IQ70"), ("all  UQW-VDB ", M_ALL, "VDB"), ("high UQW-VDB ", HIGH, "VDB"), ("all  UQW-IQ70", M_ALL, "IQ70"), ("m=1  UQW-VDB ", (1.0,), "VDB")):
            print("  ", name, st(grp(lambda n: np.mean([t2(n, "UQW", m) - t2(n, b, m) for m in ms]))))


if __name__ == "__main__":
    main()
