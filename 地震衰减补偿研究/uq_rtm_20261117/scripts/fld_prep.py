"""Stage 4 preprocessing (PLAN_STAGE4_ZH.md): Volve seq-79 / cable-5081 shot gathers and the VTI-2008 Vp
section along the line, written as raw little-endian binaries + JSON for MATLAB (cache/fld/).
usage: fld_prep.py <shots.segy> <vp.segy> <outdir>"""
import sys, os, json, math
import numpy as np

shots_fn, vel_fn, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
A = math.radians(14); X0, Y0 = 435050.02, 6478563.53
il = lambda x, y: -(x - X0) * math.cos(A) + (y - Y0) * math.sin(A)          # along 284 deg (m)
cl = lambda x, y: (x - X0) * math.sin(A) + (y - Y0) * math.cos(A)           # along 14 deg (m)
DX = 10.0; XMIN, XMAX, ZMAX = -3600.0, 2600.0, 3500.0; NPML = 25
NXI = int(round((XMAX - XMIN) / DX)) + 1; NZI = int(round(ZMAX / DX)) + 1
F = np.arange(4.0, 25.0 + 1e-9, 0.25); NT = 1001; DT = 0.004                  # 0-4 s


def ibm2f(b):
    """IBM 32-bit floats (big-endian uint32 array) -> float32."""
    b = b.astype(np.uint32); s = np.where(b >> 31, -1.0, 1.0); e = ((b >> 24) & 0x7F).astype(np.int32) - 64
    m = (b & 0x00FFFFFF).astype(np.float64) / float(1 << 24)
    return (s * m * np.power(16.0, e)).astype(np.float32)


def segy(fn):
    raw = np.memmap(fn, dtype=np.uint8, mode='r'); ns = int(np.frombuffer(raw[3220:3222].tobytes(), '>i2')[0])
    tsz = 240 + 4 * ns; n = (raw.size - 3600) // tsz
    return np.ndarray((n, tsz), np.uint8, raw[3600:3600 + n * tsz]), ns


def hdr(H, o, t='>i4'):
    w = 4 if t == '>i4' else 2
    return np.frombuffer(np.ascontiguousarray(H[:, o:o + w]).tobytes(), t).astype(np.int64)


# ---- shots ----
H, ns = segy(shots_fn)
fldr, sx, sy, gx, gy, wd = (hdr(H, o) for o in (8, 72, 76, 80, 84, 64))
scal = hdr(H, 70, '>i2'); f = np.where(scal < 0, -1.0 / scal, np.where(scal > 0, scal, 1.0))
sxs = np.array([il(a, b) for a, b in zip(sx * f, sy * f)]); gxs = np.array([il(a, b) for a, b in zip(gx * f, gy * f)])
scl = np.array([cl(a, b) for a, b in zip(sx * f, sy * f)]); gcl = np.array([cl(a, b) for a, b in zip(gx * f, gy * f)])
shot_ids = np.unique(fldr); s_pos = np.array([sxs[fldr == k][0] for k in shot_ids])
inside = (s_pos >= XMIN) & (s_pos <= XMAX); sel = shot_ids[inside][::2]
use = np.isin(fldr, sel) & (gxs >= XMIN) & (gxs <= XMAX)
idx = np.nonzero(use)[0]; print('shots', sel.size, 'traces', idx.size)
E = np.exp(-1j * 2 * np.pi * np.outer(F, np.arange(NT) * DT)).astype(np.complex64)   # nf x nt
D = np.empty((idx.size, F.size), np.complex64)
for c0 in range(0, idx.size, 2000):
    rr = idx[c0:c0 + 2000]
    tr = ibm2f(np.frombuffer(np.ascontiguousarray(H[rr, 240:240 + 4 * NT]).tobytes(), '>u4').reshape(len(rr), NT))
    D[c0:c0 + len(rr)] = tr @ E.T
shot_of = np.searchsorted(sel, fldr[idx]) + 1                                  # 1-based shot index
src_x = np.array([sxs[fldr == k][0] for k in sel]); src_ix = np.rint((src_x - XMIN) / DX).astype(int) + NPML + 1
rec_ix = np.rint((gxs[idx] - XMIN) / DX).astype(int) + NPML + 1
rec_iz = np.rint(wd[idx] / DX).astype(int) + NPML + 1                           # seabed node (water depth, m)
src_iz = np.full(sel.size, NPML + 2)                                            # 6 m -> first node below z = 0 (10 m)
D.view(np.float32).tofile(os.path.join(out, 'D.bin'))
for nm, v in dict(shot_of=shot_of, rec_ix=rec_ix, rec_iz=rec_iz, src_ix=src_ix, src_iz=src_iz).items():
    v.astype(np.int32).tofile(os.path.join(out, nm + '.bin'))
# ---- velocity: inlines 5079 and 5081 ----
V, nsv = segy(vel_fn); vil = hdr(V, 188); vxl = hdr(V, 192); cx = hdr(V, 180) * 0.01; cy = hdr(V, 184) * 0.01
secs = []
for L in (5079, 5081):
    r = np.nonzero(vil == L)[0]; r = r[np.argsort(vxl[r])]
    v = ibm2f(np.frombuffer(np.ascontiguousarray(V[r, 240:240 + 4 * nsv]).tobytes(), '>u4').reshape(len(r), nsv))
    secs.append((np.array([il(a, b) for a, b in zip(cx[r], cy[r])]), v))
s0, v0 = secs[0]; s1, v1 = secs[1]
vel = 0.5 * (v0 + np.array([np.interp(s0, s1[::-1] if s1[0] > s1[-1] else s1, (v1[::-1] if s1[0] > s1[-1] else v1)[:, k]) for k in range(nsv)]).T)
order = np.argsort(s0); s0, vel = s0[order], vel[order]
xg = XMIN + (np.arange(NXI + 2 * NPML) - NPML) * DX; zg = (np.arange(NZI + 2 * NPML) - NPML) * DX
zv = np.arange(nsv) * 10.0
vx = np.array([np.interp(xg, s0, vel[:, k]) for k in range(nsv)])              # nz_v x nx
vgrid = np.array([np.interp(np.clip(zg, 0, zv[-1]), zv, vx[:, j]) for j in range(xg.size)]).T   # nz x nx
vgrid.astype(np.float32).tofile(os.path.join(out, 'vp.bin'))
meta = dict(dx=DX, xmin=XMIN, xmax=XMAX, zmax=ZMAX, npml=NPML, nxi=NXI, nzi=NZI, nx=NXI + 2 * NPML, nz=NZI + 2 * NPML,
            f=F.tolist(), nt=NT, dt=DT, nshot=int(sel.size), ntrace=int(idx.size), fldr=sel.tolist(),
            src_cross_m=[float(np.median(scl[use])), float(scl[use].min()), float(scl[use].max())],
            rec_cross_m=[float(np.median(gcl[use])), float(gcl[use].min()), float(gcl[use].max())],
            water_depth_m=[int(wd[idx].min()), int(wd[idx].max())],
            vp_top=float(vgrid[NPML, NPML + NXI // 2]), vp_range=[float(vgrid.min()), float(vgrid.max())],
            vel_line_s=[float(s0.min()), float(s0.max())])
json.dump(meta, open(os.path.join(out, 'meta.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in meta.items() if k not in ('f', 'fldr')}, indent=1))
