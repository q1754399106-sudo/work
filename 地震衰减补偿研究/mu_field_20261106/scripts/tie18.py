"""Well ties of experiment 18, exactly as experiment 13's tie13.py, on the new mu = 0.25 sections and on experiment
13's cached comparators.

usage: tie18.py check   comparators only: reproduces experiment 13's T2 values, prints the largest difference
       tie18.py run     all methods; writes results/tie.json; prints no scores
"""
import json
import sys

from env18 import L, np, semi, lib, tie, run13, git_text, E13, CACHE13, CACHE, OUT

TW = semi.TW
COMPARATORS = {"A": ("Input", "B2_SCa"), "B": ("QSCa", "QSCa_tuned")}


def main(mode):
    logs = run13.well_logs(); G, mid = tie.groups(logs); h, x, y, il, xl = L.cube()
    cm = semi.curtain_cm(); fp = float(cm["fpeak"]); sec = {}

    def section(ln):
        if ln not in sec:
            s = {}
            for tag, names in COMPARATORS.items():
                z = np.load(CACHE13 / f"IL_{ln}_{tag}.npz")
                s.update({m: z[m] for m in names})
            if mode == "run":
                z = np.load(CACHE / f"IL_{ln}.npz"); s.update({m: z[m] for m in z.files})
            sec[ln] = s
        return sec[ln]

    rows = {}
    for n, g in logs.items():
        w, idx = tie.path(g, x, y, il, xl)
        rf = np.zeros(len(semi.TF)); tm = np.round(0.5 * (g["t"][1:] + g["t"][:-1]) - 0.0005, 6)
        i0 = int(round((tm[0] - semi.TF[0]) / L.wq.FINE)); rf[i0:i0 + len(g["rc"])] = g["rc"]
        ref = semi.ricker_ref(rf, fp); ws = semi.win(TW, g["t0"], g["t1"])
        names = list(section(int(il[idx[0, 0]]))); sc = {}
        for m in names:
            tr = np.zeros(len(TW))
            for ii, r in zip(np.where(w)[0], idx):
                tr[ii] = np.mean([section(int(il[j]))[m][ii, (xl[j] - L.XL0) // L.XSTEP] for j in r])
            sc[m] = dict(zip(("T2", "shift_ms"), lib.tie(tr, ref, TW, ws)))
        rows[n] = dict(span=[g["t0"], g["t1"]], n_traces=int(len(np.unique(idx))), scores=sc)
        print(n, "done", flush=True)
    E = json.loads(git_text(f"{E13}/results/tie.json"))
    diff = max(abs(rows[n]["scores"][m]["T2"] - E["wells"][n]["scores"][m]["T2"]) for n in rows for tag in COMPARATORS for m in COMPARATORS[tag])
    print("largest T2 difference from experiment 13 on the comparators:", diff, flush=True)
    assert G == E["groups"], "well groups differ from experiment 13"
    if mode == "run":
        OUT.mkdir(exist_ok=True)
        (OUT / "tie.json").write_text(json.dumps(dict(groups=G, wells=rows, check_vs_e13=diff), indent=1), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
