"""Repo lasso solver (benchmark.solve, unchanged: weight 1/2, tolerance as crit_kappa.py, warm start from weak compensation) on the operator of wang_op.py.

Criteria from the same solutions: J (the QScan lasso objective), E1 (eq. 17), CD (chengD), C2 (l_0.1); selection rules identical to scan_q2.m.
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; RR = HERE.parents[1]
sys.path.insert(0, str(RR / "qscan_anchor_theory_20261030" / "scripts")); sys.path.insert(0, str(HERE))
from env15 import np, sc11, bm
from scipy import sparse
import wang_op
bm.MAXITER = 20000
DT = 0.002; VGRID = np.exp(np.linspace(np.log(1 / 400), np.log(1 / 20), 21)); TCRIT = 0.25; BAND = (5.0, 90.0)


def e17(O):
    n = int(round(TCRIT / DT)); nf = int(round(1 / DT)); f = np.fft.rfftfreq(nf, DT); b = (f >= BAND[0]) & (f <= BAND[1])
    Ds = np.abs(np.fft.rfft(O[:n], nf, axis=0)).mean(axis=1); Dd = np.abs(np.fft.rfft(O[-n:], nf, axis=0)).mean(axis=1)
    return float((Ds[b].sum() - Dd[b].sum()) / max(Ds[b].sum(), 1e-30))


def cheng_d(O):
    n = O.shape[0]; m = int(round(n / 3)); h = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(m) / (m - 1)); f = np.fft.rfftfreq(256, DT)
    As = np.sqrt(np.mean(np.abs(np.fft.rfft(O[:m] * h[:, None], 256, axis=0)) ** 2, axis=1))
    Ad = np.sqrt(np.mean(np.abs(np.fft.rfft(O[-m:] * h[:, None], 256, axis=0)) ** 2, axis=1)); b = (f >= 10) & (f <= 60)
    return float(abs(np.polyfit(f[b], np.log(np.maximum(Ad[b], 1e-30) / np.maximum(As[b], 1e-30)), 1)[0]))


def gridmin(g, c):
    m = int(np.argmin(c)); v = float(g[m])
    if 0 < m < len(g) - 1:
        a = np.polyfit(g[m - 1:m + 2], c[m - 1:m + 2], 2)
        if a[0] > 0:
            v = float(np.clip(-a[1] / (2 * a[0]), g[m - 1], g[m + 1]))
    return v


def crossing(g, e):
    for k in range(len(e) - 1):
        if e[k] > 0 and e[k + 1] <= 0:
            return float(g[k] + (g[k + 1] - g[k]) * e[k] / (e[k] - e[k + 1])), True
        if k > 0 and abs(e[k]) < abs(e[k - 1]) and abs(e[k + 1]) >= abs(e[k]):
            a = np.polyfit(g[k - 1:k + 2], np.abs(e[k - 1:k + 2]), 2); v = float(g[k])
            if a[0] > 0:
                v = float(np.clip(-a[1] / (2 * a[0]), g[k - 1], g[k + 1]))
            return v, False
    return float("nan"), False


def prominent(g, J):
    best = -np.inf; v = float("nan")
    for m in range(1, len(g) - 1):
        if J[m] < J[m - 1] and J[m] <= J[m + 1]:
            p = min(J[:m].max(), J[m + 1:].max()) - J[m]
            if p > best:
                best = p; a = np.polyfit(g[m - 1:m + 2], J[m - 1:m + 2], 2); v = float(g[m])
                if a[0] > 0:
                    v = float(np.clip(-a[1] / (2 * a[0]), g[m - 1], g[m + 1]))
    return v


def scan(d, fw):
    n, nx = d.shape; live = np.ones(nx, bool); sigma = np.array([sc11.haar_sigma(d[:, j]) for j in range(nx)])
    K = sparse.csr_matrix((0, n * nx)); s2 = sigma ** 2; prev = None; rows = []; I = np.eye(n)
    for v in VGRID:
        A, W = wang_op.operator(n, DT, fw, 1 / v); scale = np.linalg.norm(A, axis=0); B = A / scale
        X, info = bm.solve(B, I, d, sigma, scale, K, 0.5, 0.0, live, start=prev); prev = X; r = X / scale[:, None]; O = W @ r
        pen1 = np.sqrt(2 * np.log(n)) * sigma[None, :] * np.median(scale) / scale[:, None]; res = d - B @ X
        J = float(np.sum((0.5 * np.sum(res ** 2, axis=0) + 0.5 * np.sum(pen1 * np.abs(X), axis=0)) / s2))
        mx = np.maximum(np.abs(r).max(axis=0), 1e-30); C2 = float(np.mean(np.sum(np.abs(r / mx) ** 0.1, axis=0)))
        rows.append((J, e17(O), cheng_d(O), C2, bool(info["converged"]), int(info["iterations"])))
    J, E, CD, C2, cv, it = map(np.array, zip(*rows)); g = VGRID
    v1, cross = crossing(g, E)
    return dict(vJ=gridmin(g, J), v1=v1, cross=cross, vd=gridmin(g, CD), v2=prominent(g, C2), J=J.tolist(), E1=E.tolist(), CD=CD.tolist(), C2=C2.tolist(),
                converged=cv.tolist(), iterations=it.tolist(), vgrid=g.tolist())
