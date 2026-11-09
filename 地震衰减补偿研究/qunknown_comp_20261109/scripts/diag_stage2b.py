"""Post-hoc diagnostic of stage 2b (written after seeing the results; descriptive only, RESULTS_STAGE2B_ZH.md):
(1) deep-window (1.0-1.75 s) power of target, noise-free data and 14 dB noise by band; (2) S at m = 1 on noise-free
and noisy data (realisation seed 1), H1-H6."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fdcheck import np, DT, load_fd, noise, score_fd, iq_t, uq_t, ricker, lib12, SEL, G_IQ, MODELS

f = np.fft.rfftfreq(451, DT); BANDS = ((5, 30), (30, 60), (60, 90), (90, 125))
for h in MODELS:
    D0, T, cm, tw, fp = load_fd(h); Q = cm["Q"]; N = noise(D0, 14, 1); w = (tw >= 1.0) & (tw <= 1.75)
    P = lambda X: np.mean(np.abs(np.fft.rfft(X[w], 451, axis=0)) ** 2, axis=1)
    pd, pn = P(D0), P(N)
    snr = ["%d-%d Hz %+.1f dB" % (a, b, 10 * np.log10(pd[(f >= a) & (f < b)].sum() / pn[(f >= a) & (f < b)].sum())) for a, b in BANDS]
    print(h, "deep-window SNR (noise-free data / 14 dB noise):", ", ".join(snr))
    for lab, D in (("noise-free", D0), ("14 dB", D0 + N)):
        r = dict(raw=score_fd(D, T, tw)[0], IQ70=score_fd(iq_t(D, cm, Q, 1, 1, G_IQ, tw), T, tw)[0],
                 IQ30=score_fd(iq_t(D, cm, Q, 1, 1, 30, tw), T, tw)[0], UQ=score_fd(uq_t(D, cm, Q, 1, 0.447, tw, 0.0), T, tw)[0],
                 VDB=score_fd(lib12.tvw(iq_t(D, cm, Q, None, 1, 0, tw), ricker(fp), SEL["B3"]["eps"], win=SEL["B3"]["win"]), T, tw)[0])
        print("   %-10s" % lab, " ".join("%s %.3f" % kv for kv in r.items()))
