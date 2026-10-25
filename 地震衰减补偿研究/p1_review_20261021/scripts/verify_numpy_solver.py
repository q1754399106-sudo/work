import json, time
from pathlib import Path
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path('/working_dir/c_d3271ae7de8a261a/extracted/attenuation_research_handoff_20261019')

def ray_absorption(t0, x, c):
    t0 = np.asarray(t0, dtype=float).ravel()
    n = len(t0)
    vp = np.asarray(c['vp'], dtype=float).ravel()
    h = np.asarray(c['h'], dtype=float).ravel()
    Q = np.asarray(c['Q'], dtype=float).ravel()
    bound = np.concatenate([[0.0], np.cumsum(2.0 * h / vp[:3]), [np.inf]])
    tau = np.maximum(0.0, np.minimum(t0[:, None], bound[None, 1:]) - bound[None, :-1])
    tau[~np.isfinite(tau)] = 0.0
    active_any = np.any(tau > 0, axis=1)
    a0 = np.sum(tau / Q[None, :], axis=1)
    a0[~active_any] = 0.0
    ax = np.zeros(n)
    ax_small = np.zeros(n)
    p = np.zeros(n)
    ax[~active_any] = x / (vp[0] * Q[0])
    if x == 0:
        return a0, ax, ax_small, p
    idx = np.where(active_any)[0]
    tt = tau[idx, :]
    act = tt > 0
    max_v = np.max(np.where(act, vp[None, :], 0.0), axis=1)
    lo = np.zeros(len(idx))
    hi = (1.0 - 1e-12) / max_v
    v2 = vp[None, :] ** 2
    for _ in range(70):
        pm = 0.5 * (lo + hi)
        pm_v = pm[:, None] * vp[None, :]
        xx = np.sum(np.where(act, tt * v2 * pm[:, None] / np.sqrt(np.maximum(1e-30, 1.0 - pm_v ** 2)), 0.0), axis=1)
        mask = xx < x
        lo = np.where(mask, pm, lo)
        hi = np.where(mask, hi, pm)
    pm = 0.5 * (lo + hi)
    p[idx] = pm
    pm_v = pm[:, None] * vp[None, :]
    ax[idx] = np.sum(np.where(act, tt / Q[None, :] / np.sqrt(np.maximum(1e-30, 1.0 - pm_v ** 2)), 0.0), axis=1) - a0[idx]
    ps = x / np.sum(tt * v2, axis=1)
    ax_small[idx] = 0.5 * (ps ** 2) * np.sum(tt * v2 / Q[None, :], axis=1)
    return a0, ax, ax_small, p

def attenuation_matrix(a, c):
    N = int(c['nfft'])
    dt = float(c['dt'])
    nt = int(c['nt'])
    fref = float(c['fref'])
    f = np.arange(0, N // 2 + 1, dtype=float)[:, None] / (N * dt)
    realmin = np.finfo(float).tiny
    b = -np.pi * f + 1j * 2.0 * f * np.log(np.maximum(f, realmin) / fref)
    b[0, 0] = 0.0
    a_row = np.asarray(a, dtype=float).reshape(1, -1)
    H = np.exp(b @ a_row)
    H[-1, :] = np.real(H[-1, :])
    Z = np.fft.irfft(H, n=N, axis=0)
    idx = (np.arange(nt)[:, None] - np.arange(nt)[None, :]) % N
    G = np.ascontiguousarray(Z[idx, np.arange(nt)[None, :]])
    return G

def ridge_operator(G, lam):
    G = np.ascontiguousarray(G)
    n = G.shape[1]
    return np.ascontiguousarray(sla.solve(G.T @ G + lam * np.eye(n), G.T, assume_a='pos'))

def discrepancy_factor(G):
    G = np.ascontiguousarray(G)
    U, s, Vt = sla.svd(G, full_matrices=False)
    return {'U': np.ascontiguousarray(U), 'z': np.ascontiguousarray(s), 'V': np.ascontiguousarray(Vt.T), 'm': G.shape[0]}

def discrepancy_apply(F, Y, grid, target):
    u = F['U'].T @ Y
    if u.ndim == 1:
        u = u[:, None]
        Y_2d = Y[:, None]
    else:
        Y_2d = Y
    perpendicular = np.maximum(0.0, np.sum(Y_2d ** 2, axis=0) - np.sum(u ** 2, axis=0))
    grid = np.asarray(grid, dtype=float).ravel()
    nr = Y_2d.shape[1]
    rr = np.zeros((len(grid), nr))
    z2 = F['z'][:, None] ** 2
    for k, g in enumerate(grid):
        rr[k, :] = np.sqrt(np.sum(((g / (z2 + g)) * u) ** 2, axis=0) + perpendicular)
    ix = np.argmin(np.abs(rr - target), axis=0)
    lam = grid[ix]
    s = F['V'] @ ((F['z'][:, None] / (z2 + lam[None, :])) * u)
    residual = rr[ix, np.arange(nr)]
    endpoint = (ix == 0) | (ix == len(grid) - 1)
    if Y.ndim == 1:
        return s.ravel(), float(lam[0]), float(residual[0]), bool(endpoint[0])
    return s, lam, residual, endpoint

def prepare_models_v2(offsets, c):
    nt = int(c['nt'])
    cache = {}
    models = []
    for x in offsets:
        xf = float(x)
        if xf in cache:
            models.append(cache[xf])
            continue
        a, b, _, p = ray_absorption(c['t'], xf, c)
        A = attenuation_matrix(a, c)
        X = attenuation_matrix(b, c)
        G = attenuation_matrix(a + b, c)
        M = ridge_operator(X, float(c['lambda_x']))
        if sla.norm(X - np.eye(nt), 'fro') / np.sqrt(nt) < 1e-10:
            M = np.eye(nt)
        L_cov = sla.cholesky(M @ M.T + float(c['cov_floor']) * np.eye(nt), lower=True)
        K = np.ascontiguousarray(sla.solve_triangular(L_cov, M, lower=True))
        F = np.ascontiguousarray(K @ (X @ A))
        fd = discrepancy_factor(G)
        fs = discrepancy_factor(F)
        m = {
            'G': G, 'K': K, 'F': F, 'A': A, 'X': X, 'M': M,
            'a0': a, 'ax': b,
            'direct': fd, 'split': fs,
            'composition_error': float(sla.norm(G - X @ A, 'fro') / sla.norm(G, 'fro')),
            'gram_direct': np.ascontiguousarray(G.T @ G),
            'gram_split': np.ascontiguousarray(F.T @ F),
            'p_time': p
        }
        cache[xf] = m
        models.append(m)
    return models

def estimate_trace_snr(d, t, noise_window=(0.02, 0.25)):
    d = np.asarray(d, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()
    normalMAD = 0.6744897501960817
    ii = (t >= noise_window[0]) & (t <= noise_window[1])
    z = d[ii]
    sigma = np.median(np.abs(z - np.median(z))) / normalMAD
    e = z - np.mean(z)
    rho1 = float((e[:-1] @ e[1:]) / max(e @ e, np.finfo(float).tiny))
    return {'sigma': float(sigma), 'samples_used': int(len(z)), 'noise_region': z, 'lag1': rho1}

def learn_structure_operator(D, c, o):
    n, nx = D.shape
    dt = float(c['dt'])
    f = np.arange(n, dtype=float) / (n * dt)
    f = np.minimum(f, 1.0 / dt - f)
    w = np.ones(n)
    p_pass = o['pilot_pass_hz']
    p_stop = o['pilot_stop_hz']
    low = f < p_pass[0]
    high = f > p_pass[1]
    w[low] = 0.5 * (1.0 - np.cos(np.pi * np.clip((f[low] - p_stop[0]) / (p_pass[0] - p_stop[0]), 0.0, 1.0)))
    w[high] = 0.5 * (1.0 + np.cos(np.pi * np.clip((f[high] - p_pass[1]) / (p_stop[1] - p_pass[1]), 0.0, 1.0)))
    pilot = np.real(np.fft.ifft(np.fft.fft(D, axis=0) * w[:, None], axis=0))
    lags = [int(x) for x in o['prediction_lags']]
    step = int(o['prediction_step_samples'])
    centers = np.unique(np.append(np.arange(0, n, step), n - 1))
    half = int(o['prediction_window_samples']) // 2
    realmin = np.finfo(float).tiny
    rows, cols, vals = [], [], []
    t_idx = np.arange(n)
    for j in range(nx):
        neighbors = [jj for jj in (j - 1, j + 1) if 0 <= jj < nx]
        if not neighbors:
            continue
        p = len(neighbors) * len(lags)
        B = np.zeros((n, p))
        ids = np.zeros((n, p), dtype=int)
        k = 0
        for jj in neighbors:
            for lag in lags:
                tt = np.clip(t_idx + lag, 0, n - 1)
                B[:, k] = pilot[tt, jj]
                ids[:, k] = tt + jj * n
                k += 1
        C = np.zeros((len(centers), p))
        gate = np.zeros(len(centers))
        for kc, c_idx in enumerate(centers):
            i0 = max(0, c_idx - half)
            i1 = min(n, c_idx + half + 1)
            Z = B[i0:i1, :]
            y = pilot[i0:i1, j]
            ridge = float(o['prediction_ridge']) * max(np.sum(Z ** 2) / p, realmin)
            coef = sla.solve(Z.T @ Z + ridge * np.eye(p), Z.T @ y, assume_a='pos')
            nrm = sla.norm(coef)
            coef = coef * min(1.0, float(o['prediction_coefficient_norm_cap']) / max(nrm, realmin))
            r2 = 1.0 - np.sum((y - Z @ coef) ** 2) / max(np.sum(y ** 2), realmin)
            C[kc, :] = coef
            gate[kc] = np.sqrt(np.clip((r2 - float(o['prediction_r2_gate'])) / (1.0 - float(o['prediction_r2_gate'])), 0.0, 1.0))
        coef_interp = np.column_stack([np.interp(t_idx, centers, C[:, k]) for k in range(p)])
        g_interp = np.interp(t_idx, centers, gate)
        g_interp = g_interp / np.sqrt(1.0 + np.sum(coef_interp ** 2, axis=1))
        rr = t_idx + j * n
        rows.append(rr)
        cols.append(rr)
        vals.append(g_interp)
        for k in range(p):
            rows.append(rr)
            cols.append(ids[:, k])
            vals.append(-g_interp * coef_interp[:, k])
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    H = sp.csc_matrix((vals, (rows, cols)), shape=(n * nx, n * nx))
    return H

def avo_projector_v2(models):
    n = models[0]['G'].shape[0]
    nx = len(models)
    ray = np.column_stack([m['p_time'] for m in models])
    realmin = np.finfo(float).tiny
    rows, cols, vals = [], [], []
    for t in range(n):
        q = (ray[t, :] ** 2)
        q = q / max(np.max(q), realmin)
        B = np.column_stack([np.ones(nx), q])
        P = B @ sla.pinv(B)
        Rt = np.eye(nx) - P
        ids = t + np.arange(nx) * n
        rr, cc = np.meshgrid(ids, ids, indexing='ij')
        rows.append(rr.ravel())
        cols.append(cc.ravel())
        vals.append(Rt.ravel())
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    return sp.csc_matrix((vals, (rows, cols)), shape=(n * nx, n * nx))

def solve_v2(ctx, scale, route):
    n, nx = ctx['n'], ctx['nx']
    lam = scale * ctx['lambda']
    mu = scale * ctx['mu_base']
    rhs = ctx['rhsD'] if route == 'Direct' else ctx['rhsS']
    u = ctx['uD'] if route == 'Direct' else ctx['uS']
    S = np.zeros((n, nx))
    N_blocks = []
    L_blocks = []
    for j in range(nx):
        m = ctx['models'][j]
        F = m['direct'] if route == 'Direct' else m['split']
        gram = m['gram_direct'] if route == 'Direct' else m['gram_split']
        S[:, j] = F['V'] @ ((F['z'] / (F['z'] ** 2 + lam[j])) * u[:, j])
        if ctx['H'].nnz > 0:
            Nj = gram + lam[j] * np.eye(n)
            N_blocks.append(Nj)
            L_blocks.append(sla.cholesky(Nj + mu * ctx['HHblocks'][j], lower=True, check_finite=False))
    if ctx['H'].nnz > 0:
        rhs_vec = rhs.ravel(order='F')
        def matvec(v):
            V = v.reshape((n, nx), order='F')
            Y = np.column_stack([N_blocks[j] @ V[:, j] for j in range(nx)])
            return Y.ravel(order='F') + mu * (ctx['HH'] @ v)
        def precond(v):
            V = v.reshape((n, nx), order='F')
            Y = np.column_stack([sla.cho_solve((L_blocks[j], True), V[:, j], check_finite=False) for j in range(nx)])
            return Y.ravel(order='F')
        linop = spla.LinearOperator((n * nx, n * nx), matvec=matvec, dtype=float)
        M_op = spla.LinearOperator((n * nx, n * nx), matvec=precond, dtype=float)
        s, flag = spla.cg(linop, rhs_vec, x0=S.ravel(order='F'), rtol=1e-7, atol=0.0, maxiter=300, M=M_op)
        normal = float(sla.norm(matvec(s) - rhs_vec) / sla.norm(rhs_vec))
        S = np.ascontiguousarray(s.reshape((n, nx), order='F'))
    else:
        normal = 0.0
        flag = 0
    rr = 0.0
    for j in range(nx):
        if route == 'Direct':
            e = ctx['models'][j]['G'] @ S[:, j] - ctx['D'][:, j]
        else:
            e = ctx['models'][j]['F'] @ S[:, j] - ctx['Yw'][:, j]
        rr += np.sum(e ** 2)
    return S, {'normal_residual': normal, 'data_residual': float(np.sqrt(rr)), 'flag': int(flag), 'solver': 'PCG' if ctx['H'].nnz > 0 else 'SVD'}

def uq_select(D, models, prior, preset, v2_opts):
    n, nx = D.shape
    gamma = float(preset['gamma']) if prior['H'].nnz > 0 else 0.0
    Yw = np.ascontiguousarray(np.column_stack([models[j]['K'] @ D[:, j] for j in range(nx)]))
    rhsD = np.ascontiguousarray(np.column_stack([models[j]['G'].T @ D[:, j] for j in range(nx)]))
    rhsS = np.ascontiguousarray(np.column_stack([models[j]['F'].T @ Yw[:, j] for j in range(nx)]))
    uD = np.ascontiguousarray(np.column_stack([models[j]['direct']['U'].T @ D[:, j] for j in range(nx)]))
    uS = np.ascontiguousarray(np.column_stack([models[j]['split']['U'].T @ Yw[:, j] for j in range(nx)]))
    HH = (prior['H'].T @ prior['H']).tocsc()
    HHblocks = [np.ascontiguousarray(HH[j*n:(j+1)*n, j*n:(j+1)*n].toarray()) for j in range(nx)]
    ctx = {
        'D': D, 'Yw': Yw, 'models': models, 'H': prior['H'], 'HH': HH, 'HHblocks': HHblocks,
        'lambda': prior['lambda'], 'mu_base': gamma * np.median(prior['lambda']),
        'rhsD': rhsD, 'rhsS': rhsS, 'uD': uD, 'uS': uS, 'n': n, 'nx': nx
    }
    target = float(preset['noise_target_multiplier']) * np.sqrt(n * np.sum(prior['sigma'] ** 2))
    grid = list(map(float, v2_opts['scale_grid']))
    cache_res = {}
    def eval_idx(idx):
        if idx not in cache_res:
            cache_res[idx] = solve_v2(ctx, grid[idx], 'Direct')
        return cache_res[idx][1]['data_residual']
    r4 = eval_idx(4)
    if r4 <= target:
        curr = 4
        while curr + 1 < len(grid) and eval_idx(curr + 1) <= target:
            curr += 1
        below = curr
        above = curr + 1 if curr + 1 < len(grid) else None
    else:
        curr = 4
        while curr - 1 >= 0 and eval_idx(curr - 1) >= target:
            curr -= 1
        above = curr
        below = curr - 1 if curr - 1 >= 0 else None
    scales = [grid[i] for i in sorted(cache_res.keys())]
    residuals = [cache_res[i][1]['data_residual'] for i in sorted(cache_res.keys())]
    candidates = [cache_res[i][0] for i in sorted(cache_res.keys())]
    records = [cache_res[i][1] for i in sorted(cache_res.keys())]
    if below is not None and above is not None and below < above:
        lo = grid[below]
        hi = grid[above]
        for _ in range(int(v2_opts['refinement_max_steps'])):
            if min(abs(r / target - 1.0) for r in residuals) <= float(v2_opts['residual_tolerance']):
                break
            sc = np.sqrt(lo * hi)
            S, r = solve_v2(ctx, sc, 'Direct')
            scales.append(sc)
            residuals.append(r['data_residual'])
            candidates.append(S)
            records.append(r)
            if r['data_residual'] < target:
                lo = sc
            else:
                hi = sc
    ix = int(np.argmin([abs(r - target) for r in residuals]))
    scale = scales[ix]
    direct = candidates[ix]
    split, sr = solve_v2(ctx, scale, 'Split')
    endpoint = int(scale == float(v2_opts['scale_grid'][0]) or scale == float(v2_opts['scale_grid'][-1]))
    mismatch = float(abs(residuals[ix] / target - 1.0))
    return {'direct': direct, 'split': split, 'scale': scale, 'target': target,
            'endpoint': endpoint, 'residual_mismatch': mismatch,
            'direct_record': records[ix], 'split_record': sr}

def make_reference(x, c):
    _, _, _, p = ray_absorption(c['event_t'], float(x), c)
    nt = int(c['nt'])
    r = np.zeros(nt)
    coeff = np.zeros(3)
    s = np.zeros(nt)
    vp_arr = np.asarray(c['vp'], dtype=float)
    vs_arr = np.asarray(c['vs'], dtype=float)
    rho_arr = np.asarray(c['rho'], dtype=float)
    t_arr = np.asarray(c['t'], dtype=float)
    ev_t = np.asarray(c['event_t'], dtype=float)
    fpeak = float(c['fpeak'])
    dt = float(c['dt'])
    for j in range(3):
        vp = np.mean(vp_arr[j:j+2])
        vs = np.mean(vs_arr[j:j+2])
        rho = np.mean(rho_arr[j:j+2])
        dvp = (vp_arr[j+1] - vp_arr[j]) / vp
        dvs = (vs_arr[j+1] - vs_arr[j]) / vs
        drho = (rho_arr[j+1] - rho_arr[j]) / rho
        A = 0.5 * (dvp + drho)
        B = 0.5 * dvp - 2.0 * ((vs / vp) ** 2) * (drho + 2.0 * dvs)
        coeff[j] = A + B * ((p[j] * vp_arr[j]) ** 2)
        z = np.pi * fpeak * (t_arr - ev_t[j])
        s += coeff[j] * (1.0 - 2.0 * (z ** 2)) * np.exp(-(z ** 2))
        r[int(round(ev_t[j] / dt))] = coeff[j]
    return s, r, coeff

def make_stability_scene(name, c):
    nx = 11
    offsets = np.arange(0, 1501, 150, dtype=float)
    lateral = np.arange(nx, dtype=float) * 50.0
    if name.endswith('Section'):
        offsets = np.full(nx, 750.0)
    n = int(c['nt'])
    truth = np.zeros((n, nx))
    events, coeff, templates, ang = [], [], [], []
    t_arr = np.asarray(c['t'], dtype=float)
    ev_t = np.asarray(c['event_t'], dtype=float)
    fpeak = float(c['fpeak'])
    vp_arr = np.asarray(c['vp'], dtype=float)
    for j in range(nx):
        _, _, a = make_reference(offsets[j], c)
        a = a.copy()
        times = ev_t.copy()
        if name == 'FaultSection':
            times = times + 0.008 * np.sin(j * np.pi / 10.0)
            if j >= 6:
                times[1:3] = times[1:3] + 0.04
        elif name == 'PinchoutSection':
            times = np.array([ev_t[0], ev_t[1], ev_t[1] + 0.030 * max(0.0, (8 - j) / 8.0)])
            a = np.array([a[0], a[1], -a[1]])
        elif name == 'WeakSection':
            times[1] = times[1] + 0.012 * np.sin(j * np.pi / 10.0)
            a[1] = 0.08 * a[1]
        W = np.zeros((n, 3))
        for k in range(3):
            z = np.pi * fpeak * (t_arr - times[k])
            W[:, k] = (1.0 - 2.0 * (z ** 2)) * np.exp(-(z ** 2))
        truth[:, j] = W @ a
        events.append(times)
        coeff.append(a)
        templates.append(W)
        _, _, _, p = ray_absorption(ev_t, offsets[j], c)
        ang.append((p * vp_arr[:3]) ** 2)
    return {'name': name, 'truth': np.ascontiguousarray(truth), 'offsets': offsets,
            'lateral': lateral, 'events': events, 'coeff': coeff, 'templates': templates, 'sin2': ang}

def score_stability(S, scene, c):
    truth = scene['truth']
    t_arr = np.asarray(c['t'], dtype=float)
    ev_t = np.asarray(c['event_t'], dtype=float)
    met = {}
    met['RelativeError'] = float(sla.norm(S - truth, 'fro') / sla.norm(truth, 'fro'))
    ii = (t_arr >= 0.28) & (t_arr <= 1.08)
    met['SignalWindowError'] = float(sla.norm(S[ii, :] - truth[ii, :], 'fro') / sla.norm(truth[ii, :], 'fro'))
    nn = (t_arr >= 0.02) & (t_arr <= 0.25)
    met['NoiseWindowRMS'] = float(np.sqrt(np.mean(S[nn, :] ** 2)))
    nx = S.shape[1]
    amp = np.zeros((3, nx))
    amptrue = np.zeros((3, nx))
    theta = np.zeros((3, nx))
    unstable = np.zeros((3, nx), dtype=bool)
    for j in range(nx):
        W = scene['templates'][j]
        amp[:, j] = sla.pinv(W) @ S[:, j]
        amptrue[:, j] = scene['coeff'][j]
        theta[:, j] = scene['sin2'][j]
        if np.linalg.cond(W.T @ W) > 1e6:
            unstable[:, j] = True
    valid = ~unstable
    met['EventAmplitudeError'] = float(sla.norm(amp[valid] - amptrue[valid]) / max(sla.norm(amptrue[valid]), np.finfo(float).tiny))
    met['UnidentifiableAmplitudes'] = int(np.sum(unstable))
    met['AVOInterceptError'] = np.nan
    met['AVOGradientError'] = np.nan
    if scene['name'].endswith('AVO'):
        est = np.zeros((2, 3))
        ref = np.zeros((2, 3))
        for k in range(3):
            B = np.column_stack([np.ones(nx), theta[k, :]])
            est[:, k] = sla.lstsq(B, amp[k, :])[0]
            ref[:, k] = sla.lstsq(B, amptrue[k, :])[0]
        met['AVOInterceptError'] = float(sla.norm(est[0, :] - ref[0, :]) / sla.norm(ref[0, :]))
        met['AVOGradientError'] = float(sla.norm(est[1, :] - ref[1, :]) / sla.norm(ref[1, :]))
    met['WeakEventError'] = np.nan
    met['FaultNeighborhoodError'] = np.nan
    met['PinchoutTailError'] = np.nan
    if scene['name'] == 'WeakSection':
        met['WeakEventError'] = float(sla.norm(amp[1, :] - amptrue[1, :]) / sla.norm(amptrue[1, :]))
    elif scene['name'] == 'FaultSection':
        k = slice(4, 8)
        met['FaultNeighborhoodError'] = float(sla.norm(S[:, k] - truth[:, k], 'fro') / sla.norm(truth[:, k], 'fro'))
    elif scene['name'] == 'PinchoutSection':
        k = slice(8, 11)
        t_mask = (t_arr > ev_t[1] - 0.05) & (t_arr < ev_t[1] + 0.06)
        met['PinchoutTailError'] = float(sla.norm(S[t_mask, k] - truth[t_mask, k], 'fro') / sla.norm(truth, 'fro'))
    return met
