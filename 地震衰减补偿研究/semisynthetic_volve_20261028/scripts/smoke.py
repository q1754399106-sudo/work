"""Smoke test before the plan is frozen: every method and score on a 31-trace crop with Gaussian noise of
a seed that no dataset uses. Prints only run status and timing, not scores."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run, semi  # noqa: E402
from semi import np, TW  # noqa: E402

cm = semi.curtain_cm(); B = np.load(run.BUILD); sl = slice(85, 116); semi.CENTRE = 15
E = semi.gaussian_noise(B["E0"][:, sl], np.random.default_rng(999))
D = B["sig_c1.0"][:, sl] + E
t0 = time.monotonic(); out, info = run.methods(D, cm, 1.0)
x = np.asarray(B["src_c1.0"])[:, sl]
sc = {m: semi.scores(A, x, x, x[:, 0]) for m, A in out.items()}
assert all(np.isfinite(v) for s in sc.values() for v in s.values()), "non-finite score"
print("smoke ok", sorted(out), f"{time.monotonic() - t0:.0f}s", {k: v.get("converged") for k, v in info.items() if isinstance(v, dict)})
