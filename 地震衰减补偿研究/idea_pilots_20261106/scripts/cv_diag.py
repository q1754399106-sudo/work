"""One diagnostic of the split criterion: why does it prefer high gain? Tests on 4 development datasets:
(1) real noise halves vs phase-randomised (truly independent) halves; (2) full band vs lowpass at f_true in the criterion;
(3) the target of the risk: band-limited ideal vs full-band src. Exploratory."""
import json, sys
sys.path.insert(0, r"D:\ChatGPT\研究\地震衰减补偿研究\mu_field_20261106\scripts")
from env18 import L, np, semi, git_text, RES

lib12 = L.lib12
META = json.loads(git_text("semisynth_fault_q_20261030/analysis/datasets.json"))
DATA = RES / "semisynth_fault_q_20261030" / "cache" / "data"
GRID = [30, 40, 50, 60, 70, 80, 90, 100]; REF = 120
cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); rng = np.random.default_rng(0)


def load(d):
    return {k: v.astype(float) for k, v in np.load(DATA / f"{d}.npz").items()}


def randomise(n):
    F = np.fft.rfft(n, axis=0); ph = np.exp(2j * np.pi * rng.random(F.shape)); ph[0] = 1
    return np.fft.irfft(F * ph, n.shape[0], axis=0)


for d, o in (("ZA_c1.05_r0", "ZA_c1.05_r1"), ("YA_c1.05_r0", "YA_c1.05_r1"), ("YA_c1.4_r0", "YA_c1.4_r1"), ("YA_c0.75_r0", "YA_c0.75_r1")):
    a, b = load(d), load(o); c = META[d]["c_eff_cum"]; ft = META[d]["f_true"]; x, n1 = a["signal"], a["noise"]
    n2 = b["noise"] * np.sqrt(np.mean(n1 ** 2) / np.mean(b["noise"] ** 2))
    iq = lambda D, g: lib12.iq(D, cm, Q, c, c, g); lp = lambda X: semi.lowpass(X, ft)
    F1, F2 = np.fft.rfft(n1, axis=0), np.fft.rfft(n2, axis=0)
    coh = np.abs(np.mean(F1 * np.conj(F2), axis=1)) / np.sqrt(np.mean(np.abs(F1) ** 2, axis=1) * np.mean(np.abs(F2) ** 2, axis=1)); f = np.fft.rfftfreq(n1.shape[0], semi.DT)
    print(d, f"f_true {ft:.0f}; noise-line coherence median 10-60 Hz {np.median(coh[(f > 10) & (f < 60)]):.2f}, 60-120 Hz {np.median(coh[(f > 60) & (f < 120)]):.2f}; expected if independent about {1 / np.sqrt(n1.shape[1]):.2f}", flush=True)
    full = x + (n1 + n2) / np.sqrt(2)
    S = [lib12.score(iq(full, g), a["ideal"])["S"] for g in GRID]
    mi = [float(np.mean((iq(full, g) - a["ideal"]) ** 2)) for g in GRID]; ms = [float(np.mean((iq(full, g) - a["src"]) ** 2)) for g in GRID]
    print("  S", np.round(S, 3).tolist(), "best S", GRID[int(np.argmax(S))], "| best MSE ideal", GRID[int(np.argmin(mi))], "| best MSE src", GRID[int(np.argmin(ms))])
    for name, (m1, m2) in dict(real=(n1, n2), randomised=(randomise(n1), randomise(n2))).items():
        y1, y2 = x + np.sqrt(2) * m1, x + np.sqrt(2) * m2; r1, r2 = iq(y1, REF), iq(y2, REF); cf, cl = [], []
        for g in GRID:
            u, v = iq(y1, g), iq(y2, g)
            cf.append(float(np.mean((u - r2) ** 2) + np.mean((v - r1) ** 2)))
            cl.append(float(np.mean((lp(u) - lp(r2)) ** 2) + np.mean((lp(v) - lp(r1)) ** 2)))
        print(f"  {name}: pick full band {GRID[int(np.argmin(cf))]}, pick lowpass at f_true {GRID[int(np.argmin(cl))]}; S at picks {S[int(np.argmin(cf))]:.3f}, {S[int(np.argmin(cl))]:.3f}", flush=True)
