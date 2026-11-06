"""Exploratory pilot, used development data only (structure A). Can an image-domain split criterion pick the inverse-Q
gain without truth?  C(g) = || iq_g(half 1) - iq_ref(half 2) ||^2: half 2's nearly unregularised inverse is a noisy,
nearly unbiased target with noise independent of half 1 (Noise2Noise logic), so C(g) should track the true risk.
A: two simulated half stacks with independent real noise. B: adjacent traces of the post-stack section as the split.
True Q (c_eff_cum) is used throughout, so only the gain is being selected."""
import json, sys
sys.path.insert(0, r"D:\ChatGPT\研究\地震衰减补偿研究\mu_field_20261106\scripts")
from env18 import L, np, semi, git_text, RES

lib12 = L.lib12; TW = semi.TW
META = json.loads(git_text("semisynth_fault_q_20261030/analysis/datasets.json"))
DATA = RES / "semisynth_fault_q_20261030" / "cache" / "data"
GRID = [20, 30, 40, 50, 60, 70, 80, 90, 100]; REF = 120
cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float)


def load(d):
    return {k: v.astype(float) for k, v in np.load(DATA / f"{d}.npz").items()}


rows = []
for var in "ZY":
    for c0 in (0.75, 1.05, 1.4):
        for r in (0, 1):
            d = f"{var}A_c{c0}_r{r}"; a = load(d); b = load(f"{var}A_c{c0}_r{1 - r}"); c = META[d]["c_eff_cum"]
            x, n1 = a["signal"], a["noise"]; n2 = b["noise"] * np.sqrt(np.mean(n1 ** 2) / np.mean(b["noise"] ** 2))
            full = x + (n1 + n2) / np.sqrt(2); y1, y2 = x + np.sqrt(2) * n1, x + np.sqrt(2) * n2
            iq = lambda D, g: lib12.iq(D, cm, Q, c, c, g)
            r1, r2 = iq(y1, REF), iq(y2, REF)
            D = a["D"]; lags = semi.local_lags(D); rD = iq(D, REF); n, nx = D.shape; idx = np.arange(n, dtype=float)
            rb = np.stack([semi.sinc_shift(rD[:, j + 1], idx + lags[:, j]) for j in range(nx - 1)], axis=1)
            sF, sD, cA, cB, risk = [], [], [], [], []
            for g in GRID:
                o = iq(full, g); sF.append(lib12.score(o, a["ideal"])["S"]); risk.append(float(np.mean((o - a["ideal"]) ** 2)))
                u = iq(D, g); sD.append(lib12.score(u, a["ideal"])["S"])
                cA.append(float(np.mean((iq(y1, g) - r2) ** 2) + np.mean((iq(y2, g) - r1) ** 2)))
                cB.append(float(np.mean((u[:, :-1] - rb) ** 2)))
            rows.append(dict(id=d, sF=sF, sD=sD, cA=cA, cB=cB, risk=risk))
            print(d, "S", [round(v, 3) for v in sF], "| A pick", GRID[int(np.argmin(cA))], "B pick", GRID[int(np.argmin(cB))],
                  "oracle S", GRID[int(np.argmax(sF))], "oracle MSE", GRID[int(np.argmin(risk))], flush=True)
G = np.array(GRID)
print("\nmean S on the full data per gain:", dict(zip(GRID, np.round(np.mean([r["sF"] for r in rows], axis=0), 3).tolist())))
for name, crit, s in (("A half stacks", "cA", "sF"), ("B adjacent traces", "cB", "sD")):
    p = [int(np.argmin(r[crit])) for r in rows]; val = [r[s][i] for r, i in zip(rows, p)]; orc = [max(r[s]) for r in rows]
    fix = {g: np.mean([r[s][GRID.index(g)] for r in rows]) for g in (30, 70)}
    print(f"{name}: picks {G[p].tolist()}; mean S picked {np.mean(val):.3f}, oracle per dataset {np.mean(orc):.3f}, fixed 30 dB {fix[30]:.3f}, fixed 70 dB {fix[70]:.3f}; picked - fixed 70 per dataset {np.round(np.array(val) - np.array([r[s][GRID.index(70)] for r in rows]), 3).tolist()}")
    rho = [np.corrcoef(r[crit], r["risk"])[0, 1] for r in rows]; print("   correlation of criterion with true MSE per dataset (median):", round(float(np.median(rho)), 3))
print("oracle-gain spread across datasets (S):", G[[int(np.argmax(r["sF"])) for r in rows]].tolist())
json.dump(dict(grid=GRID, rows=rows), open(r"C:\Users\qq175\AppData\Local\Temp\claude\cb810a85168e9b1c\0777d816-dbf9-4ae1-b342-2f76f5565fae\scratchpad\cv_pilot.json", "w"))
