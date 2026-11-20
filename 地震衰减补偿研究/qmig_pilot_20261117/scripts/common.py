"""Shared constants, model and Helmholtz solver for the Q-compensated migration pilot (exploratory, synthetic).

Convention: p(t) = sum_f P(f) exp(-i w t); outgoing waves exp(+i k r); Im k > 0 attenuates.
Complex slowness of an operator (d, kappa):  (1/c0) [1 - d qi ln(f/F0)/pi] [1 + i kappa qi/2],  qi = 1/Q.
d = 1, kappa = 1: constant-Q (Kolsky-Futterman) medium; d = kappa = 0: acoustic at the reference frequency F0;
kappa < 0: amplitude-boosting (compensating) operator with the same dispersion.
"""
import os, sys
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get("QMIG_CACHE", ROOT / "cache")); CACHE.mkdir(parents=True, exist_ok=True)

H = 8.0; NX, NZ = 301, 176; NPML = 25
F0 = 20.0; FPEAK = 18.0
TREC = 1.6; DF = 1.0 / TREC; JF = np.arange(5, 69); FREQ = JF * DF            # 3.125 ... 42.5 Hz, 64 frequencies
ZSR = 2                                                                       # source / receiver depth index (16 m)
IXR = np.arange(0, NX, 2)                                                     # 151 receivers, 16 m apart
IXS = IXR[6:150:6]                                                            # 24 shots, 96 m apart
MS = (0.5, 0.7, 1.0, 1.4, 2.0)                                                # Q multipliers of the migration model
Q_BG, Q_LENS = 80.0, 18.0
SEED = 20261117


def ricker(f):
    return 2 / np.sqrt(np.pi) * f ** 2 / FPEAK ** 3 * np.exp(-(f / FPEAK) ** 2)


def model():
    """True c0 (m/s at F0) and 1/Q on the (NX, NZ) grid; smooth migration versions."""
    rng = np.random.default_rng(SEED)
    x = np.arange(NX) * H; z = np.arange(NZ) * H; X, Z = np.meshgrid(x, z, indexing="ij")
    c = 2000 + 0.8 * Z
    hgt = 60 * np.exp(-((X - 1200) / 500) ** 2) - 40 * (X > 1500)            # anticline and a fault (40 m throw)
    zs = Z + hgt                                                              # structural depth
    edges = [170.0]
    while edges[-1] < 1700:
        edges.append(edges[-1] + rng.uniform(24, 72))
    pert = rng.uniform(-0.08, 0.08, len(edges)); pert[0] = 0
    idx = np.clip(np.searchsorted(edges, zs, side="right"), 0, len(edges)) - 1
    c = c * (1 + np.where(idx >= 0, pert[np.maximum(idx, 0)], 0.0))
    e = np.sqrt(((X - 1200) / 420) ** 2 + ((Z - 450) / 110) ** 2)            # gas lens
    wl = np.clip((1.15 - e) / 0.3, 0, 1); wl = wl * wl * (3 - 2 * wl)
    c = c * (1 - 0.06 * wl)
    qi = (1 - wl) / Q_BG + wl / Q_LENS
    sm = lambda a: gaussian_filter(a, 5, mode="nearest")
    cs = 1 / sm(1 / c); qis = sm(qi)
    return dict(c=c, qi=qi, cs=cs, qis=qis, x=x, z=z)


def _d1(n, s_int, s_half):
    cf = np.array([1, -27, 27, -1]) / (24 * H)
    Df = sp.diags(cf, [-1, 0, 1, 2], shape=(n, n)); Db = sp.diags(cf, [-2, -1, 0, 1], shape=(n, n))
    return sp.diags(1 / s_int) @ Db @ sp.diags(1 / s_half) @ Df


def _pad(a):
    return np.pad(a, NPML, mode="edge")


def laplacian(f, cpml=2500.0):
    w = 2 * np.pi * f; L = NPML * H; smax = -3 * cpml * np.log(1e-4) / (2 * L)
    out = []
    for n in (NX, NZ):
        nt = n + 2 * NPML; i = np.arange(nt, dtype=float)
        dist = lambda p: np.maximum(np.maximum(NPML - p, p - (nt - 1 - NPML)), 0) * H
        s = lambda p: 1 + 1j * smax * (dist(p) / L) ** 2 / w
        out.append(_d1(nt, s(i), s(i + 0.5)))
    nxt, nzt = NX + 2 * NPML, NZ + 2 * NPML
    return (sp.kron(out[0], sp.eye(nzt)) + sp.kron(sp.eye(nxt), out[1])).tocsc()


def slowness(c0, qi, f, d, kappa):
    return (1 / c0) * (1 - d * qi * np.log(f / F0) / np.pi) * (1 + 0.5j * kappa * qi)


def factor(lap, c0, qi, f, d, kappa):
    k2 = (2 * np.pi * f * _pad(slowness(c0, qi, f, d, kappa))) ** 2
    return spl.splu((lap + sp.diags(k2.ravel())).tocsc())


def gidx(ix, iz):
    return (np.asarray(ix) + NPML) * (NZ + 2 * NPML) + (np.asarray(iz) + NPML)


def interior(v):
    """(ntot, nrhs) -> (NX, NZ, nrhs)"""
    return v.reshape(NX + 2 * NPML, NZ + 2 * NPML, -1)[NPML:-NPML, NPML:-NPML]
