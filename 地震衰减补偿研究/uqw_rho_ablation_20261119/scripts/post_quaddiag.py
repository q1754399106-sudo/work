"""Post-hoc check (not in PLAN_ZH.md): where do the V0d and V0 outputs differ in time? Semi ZA_c0.75_r0, m 0.5 and 2.
Explains why the smoke output difference (0.26 / 0.098) gives |S(V0d) - S(V0)| ~ 1e-13: the S window is 2.30-2.80 s."""
from ablate import *

cm = semi.curtain_cm(); Q = np.asarray(cm["Q"], float); arr, meta = load("ZA_c0.75_r0"); D = arr["D"]; inv = snr_rows(D, TW, F)
for m in (0.5, 2.0):
    c = m * meta["c_eff_cum"]
    ref = U.uqw(D, cm, Q, c, SIG, TW, 1.4, 0.0, inv_snr=inv); x = uqw_var(D, cm, Q, c, SIG, TW, 1.4, inv, gl_rule(E_LO, E_HI, N_FULL), "scaled")
    d = np.abs(x - ref).max(axis=1); tol = 1e-6 * np.abs(ref).max(); w = np.where(d > tol)[0]
    print("m %.1f  output max rel %.3f; samples with diff > 1e-6 of max: t %.3f-%.3f s; max rel diff for t >= 1.4 s %.1e; in 2.30-2.80 s %.1e"
          % (m, d.max() / np.abs(ref).max(), TW[w.min()], TW[w.max()], d[TW >= 1.4].max() / np.abs(ref).max(),
             d[(TW >= 2.3) & (TW <= 2.8)].max() / np.abs(ref).max()))
