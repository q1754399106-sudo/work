"""Model-only quantities per profile, computed before any test data (PLAN_ZH.md table 2): r, uniform U, a0/a(2.10 s)."""
from core import np, lib, semi, bm, PROFILES, FAM_A
cm = semi.curtain_cm(); t = bm.vec(cm["t"]).astype(float)
print("| 剖面 | 族 | Q（4 层） | r | 均匀 U | a₀/a(2.10 s) | 差 |")
print("|---|---|---|---|---|---|---|")
for p, Q in PROFILES.items():
    c = dict(cm, Q=np.asarray(Q, float)); a = lib.absorption(c, t); a0 = float(lib.absorption(c, [1.4])[0])
    A = lambda x: float(lib.absorption(c, [x])[0]); r = (A(3.0) - A(2.1)) / (A(2.1) - A(1.2))
    U = a0 * a.mean() / (a ** 2).mean(); m = a0 / A(2.10)
    print(f"| {p} | {'A' if p in FAM_A else 'B'} | {'/'.join(f'{q:g}' for q in Q)} | {r:.2f} | {U:.3f} | {m:.3f} | {m - U:+.3f} |")
