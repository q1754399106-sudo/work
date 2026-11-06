"""Baselines on offset gathers. SC: Ma et al. (2023) eq. 17 as implemented in novelty_evidence_20261025/scripts/benchmark.py, with its frozen selection
(a = 0.125, mu_relative = 10; SC0: a = 0.5, mu = 0). Changes: per-offset dictionary (the original has one dictionary for a common-offset section); the predictor is
written with array operations instead of scipy.sparse, same stencil and same least-squares fit; the gather's offset axis plays the role of the lateral axis.
Data exchange with MATLAB by raw float64 files, column-major."""
import json, time
from pathlib import Path
import numpy as np
from gate1 import C, NX, n, DT, T, a0_ax, atten, ricker

EX = Path(__file__).resolve().parents[1] / "ex"; EX.mkdir(exist_ok=True); MAXITER = 5000; TOL = 0.02


def dictionaries(Q, fp):
    a0, ax = a0_ax(Q); z = np.pi * fp * (T[:, None] - T[None, :]); W = (1 - 2 * z ** 2) * np.exp(-z ** 2); out = []
    for i in range(NX):
        B = atten(a0 + ax[:, i]) @ W; s = np.linalg.norm(B, axis=0); B = B / s; out.append((B, B.T @ B, s))
    return W, out


class Pred:
    """K x = -x(t, j) + sum b_ij x(t + i, j + j'), same stencil as benchmark.predictor. Fitted on the data."""
    OFF = [(i, j) for j in (-2, -1, 1, 2) for i in range(-2, 3)]

    def __init__(self, D):
        nn, nx = D.shape; self.tt, self.jj = np.meshgrid(np.arange(2, nn - 2), np.arange(2, nx - 2), indexing="ij"); self.tt = self.tt.ravel(); self.jj = self.jj.ravel(); self.shape = D.shape
        X = np.column_stack([D[self.tt + i, self.jj + j] for i, j in self.OFF]); self.b = np.linalg.lstsq(X, D[self.tt, self.jj], rcond=None)[0]

    def K(self, x):
        r = -x[self.tt, self.jj].copy()
        for c, (i, j) in zip(self.b, self.OFF):
            r += c * x[self.tt + i, self.jj + j]
        return r

    def Kt(self, r):
        x = np.zeros(self.shape); x[self.tt, self.jj] -= r      # each index pair below is unique, so plain += is the adjoint
        for c, (i, j) in zip(self.b, self.OFF):
            x[self.tt + i, self.jj + j] += c * r
        return x


def solve_sc(D, sig, W, dicts, a, mu_rel):
    """D (n, NX), sig (NX). FISTA on 0.5|B x - D|^2 + pen |x| + 0.5 mu |K P x|^2, per-offset B."""
    nn, nx = D.shape; gram = [d[1] for d in dicts]; rhs = np.stack([dicts[i][0].T @ D[:, i] for i in range(nx)], 1); scale = np.stack([d[2] for d in dicts], 1)
    Ld = max(float(np.linalg.eigvalsh(g)[-1]) for g in gram); P = W
    if mu_rel > 0:
        pr = Pred(D)
        def sh(x):
            Px = P @ (x / scale); return (P.T @ pr.Kt(pr.K(Px))) / scale
        rng = np.random.default_rng(27002); v = rng.normal(size=D.shape); v /= np.linalg.norm(v)
        for _ in range(35):
            w = sh(v); nv = np.linalg.norm(w)
            if nv == 0:
                break
            v = w / nv
        Ls = float(np.sum(v * sh(v))); mu = mu_rel * Ld / max(Ls, 1e-30)
    else:
        sh = None; Ls = 0.0; mu = 0.0
    L = 1.05 * (Ld + mu * Ls); med = np.median(scale); pen = np.maximum(a * np.sqrt(2 * np.log(nn)) * sig[None, :] * med / scale, 1e-14)
    def h(z):
        o = np.stack([gram[i] @ z[:, i] for i in range(nx)], 1)
        return o + mu * sh(z) if mu else o
    x = np.zeros_like(D); y = x.copy(); th = 1.0; kkt = np.inf; it = 0
    for it in range(1, MAXITER + 1):
        grad = h(y) - rhs
        while True:
            z = y - grad / L; xx = np.sign(z) * np.maximum(np.abs(z) - pen / L, 0); dl = xx - y
            if np.sum(dl * h(dl)) > L * np.sum(dl * dl) * (1 + 1e-10):
                L *= 2; continue
            break
        nt_ = (1 + np.sqrt(1 + 4 * th * th)) / 2; yy = xx + (th - 1) / nt_ * (xx - x)
        if np.sum((y - xx) * (xx - x)) > 0:
            yy = xx.copy(); nt_ = 1.0
        x, y, th = xx, yy, nt_
        if it % 25 == 0:
            g = h(x) - rhs; rr = np.maximum(np.abs(g) - pen, 0); act = x != 0; rr[act] = np.abs(g[act] + (pen * np.sign(x))[act]); kkt = float(np.max(rr / pen))
            if kkt <= TOL:
                break
    return W @ (x / scale), dict(it=it, kkt=kkt, conv=bool(kkt <= TOL))


def write_cell(name, D, Q, fp, sig, off=C.offset):
    (EX / f"{name}_D.bin").write_bytes(np.asfortranarray(D).tobytes(order="F"))
    (EX / f"{name}_meta.json").write_text(json.dumps(dict(name=name, Q=list(map(float, Q)), fp=float(fp), sig=np.asarray(sig).tolist(), nc=int(D.shape[2]), offsets=list(map(float, off)))), encoding="utf-8")


def read_out(name, kind, nc):
    return np.fromfile(EX / f"{name}_{kind}.bin", dtype=np.float64).reshape((n, NX, nc), order="F")
