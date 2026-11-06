"""Development only (seed 901): checks the operator port against MATLAB, then one lasso scan per condition at Q 70."""
import sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from lasso_scan import np, scan, wang_op, DT
from scipy.io import loadmat
C = HERE.parent / "cache"; n = 1001
A = loadmat(C / "dev_A_fw24_Q70.mat")["A"]; Ap, _ = wang_op.operator(n, DT, 24, 70.0)
print("operator port: rel. difference", np.linalg.norm(A - Ap) / np.linalg.norm(A), flush=True)
m = loadmat(C / "dev_d_Q70_s901.mat"); d = m["dd"]; fp = float(m["fp"].ravel()[0]); ts = 1.0
for name, fw in (("known", 30.0), ("W3_0", fp), ("W3_1", 1 / np.sqrt(1 / fp**2 - np.pi * ts * 0.01 / (2 * fp)))):
    t0 = time.monotonic(); S = scan(d, fw)
    print(name, "fw %.1f Q: J %.0f e17 %.0f (cross %s) chengD %.0f C2 %.0f | conv %.2f it_med %d | %.0fs" % (fw, 1 / S["vJ"], 1 / S["v1"], S["cross"], 1 / S["vd"], 1 / S["v2"], np.mean(S["converged"]), np.median(S["iterations"]), time.monotonic() - t0), flush=True)
    print("  E1", np.round(S["E1"][::2], 2).tolist(), flush=True)
