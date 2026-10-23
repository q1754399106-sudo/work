"""
Core P1 numerical solver and calibration functions (SciPy numerical twin).
Enforces public input whitelist, zero truth access, genuine Q/fpeak posterior uncertainty propagation,
and exact normal residual verification across all branches (including NoStructure).
"""
import math, sys
from pathlib import Path
import numpy as np
import scipy.linalg as sla
import scipy.optimize as sopt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_numpy_solver as v

def smooth_box(x, w=19):
    k = np.ones(w, dtype=float) / w
    if x.ndim == 1:
        return np.convolve(x, k, mode='same')
    return np.column_stack([np.convolve(x[:, j], k, mode='same') for j in range(x.shape[1])])

def enrich_calibration_uncertainty(est):
    """Compute posterior log-(Q, f) std, correlation, and combined uncertainty from profile_losses surface."""
    qgrid = np.asarray(est['q_grid'], dtype=float)
    fgrid = np.asarray(est['frequency_grid'], dtype=float)
    losses = np.asarray(est['profile_losses'], dtype=float)
    best = float(np.min(losses))
    threshold = float(est['near_optimal_threshold'])
    temp = max(float(est["near_optimal_threshold"]) - float(est["minimum_loss"]), 0.15 * best, 1.0)
    W = np.exp(-np.maximum(0.0, losses - best) / temp)
    W /= max(float(np.sum(W)), np.finfo(float).tiny)
    FF, QQ = np.meshgrid(np.log(fgrid), np.log(qgrid))
    mq = float(np.sum(W * QQ))
    mf = float(np.sum(W * FF))
    vq = float(np.sum(W * ((QQ - mq) ** 2)))
    vf = float(np.sum(W * ((FF - mf) ** 2)))
    cov_qf = float(np.sum(W * (QQ - mq) * (FF - mf)))
    log_q_std = float(math.sqrt(max(vq, 1e-4)))
    log_f_std = float(math.sqrt(max(vf, 1e-4)))
    log_qf_corr = float(cov_qf / max(math.sqrt(vq * vf), np.finfo(float).tiny))
    q_range = est['q_sensitivity_range']
    f_range = est['frequency_sensitivity_range']
    q_span = 0.5 * math.log(max(q_range) / max(min(q_range), 1.0))
    f_span = 0.5 * math.log(max(f_range) / max(min(f_range), 1.0))
    q_unc = float(math.sqrt(log_q_std ** 2 + 0.25 * (q_span ** 2)))
    f_unc = float(math.sqrt(log_f_std ** 2 + 0.25 * (f_span ** 2)))
    est['log_q_std'] = log_q_std
    est['log_f_std'] = log_f_std
    est['log_qf_corr'] = log_qf_corr
    est['q_uncertainty'] = q_unc
    est['f_uncertainty'] = f_unc
    return est

def build_inverse_config(public, est):
    """Mirror uq_inverse_config.m: propagate calibrated Q, fpeak, and their uncertainty metrics into cm."""
    assert 'Q' not in public and 'fpeak' not in public and 'rho' not in public and 'vs' not in public
    cm = dict(public)
    cm['Q'] = [float(est['effective_Q'])] * 4
    cm['fpeak'] = float(est['estimated_frequency'])
    cm['event_t'] = list(np.cumsum(2.0 * np.asarray(public['h'], dtype=float) / np.asarray(public['vp'][:3], dtype=float)))
    cm['log_q_std'] = float(est.get('log_q_std', 0.22))
    cm['log_f_std'] = float(est.get('log_f_std', 0.08))
    cm['log_qf_corr'] = float(est.get('log_qf_corr', 0.35))
    cm['q_sensitivity_range'] = list(est.get('q_sensitivity_range', [cm['Q'][0] * 0.85, cm['Q'][0] * 1.15]))
    cm['frequency_sensitivity_range'] = list(est.get('frequency_sensitivity_range', [cm['fpeak'] * 0.92, cm['fpeak'] * 1.08]))
    q_span = 0.5 * math.log(max(cm['q_sensitivity_range']) / max(min(cm['q_sensitivity_range']), 1.0))
    f_span = 0.5 * math.log(max(cm['frequency_sensitivity_range']) / max(min(cm['frequency_sensitivity_range']), 1.0))
    cm['q_uncertainty'] = float(est.get('q_uncertainty', math.sqrt(cm['log_q_std'] ** 2 + 0.25 * (q_span ** 2))))
    cm['f_uncertainty'] = float(est.get('f_uncertainty', math.sqrt(cm['log_f_std'] ** 2 + 0.25 * (f_span ** 2))))
    return cm

def prepare_models_with_fmode(offsets, cm):
    """Build forward/inverse operators for given cm and attach singular mode centroid frequencies f_mode."""
    invmodels = v.prepare_models_v2(offsets, cm)
    nt = int(cm['nt'])
    dt = float(cm['dt'])
    freqs = np.fft.rfftfreq(nt, dt)
    for m in invmodels:
        for rk in ('direct', 'split'):
            V = m[rk]['V']
            Vf = np.abs(np.fft.rfft(V, axis=0)) ** 2
            m[rk]['f_mode'] = (freqs @ Vf) / np.maximum(np.sum(Vf, axis=0), np.finfo(float).tiny)
    return invmodels

def calibrate_q_from_gather(D, offsets, public, uq_spec):
    assert 'Q' not in public and 'fpeak' not in public and 'rho' not in public and 'vs' not in public
    traces = [idx - 1 for idx in uq_spec['calibration_trace_indices']]
    t_arr = np.asarray(public['t'], dtype=float)
    times = np.cumsum(2.0 * np.asarray(public['h'], dtype=float) / np.asarray(public['vp'][:3], dtype=float))
    nt = int(public['nt'])
    mask = np.zeros(nt, dtype=bool)
    for t0 in times:
        mask |= np.abs(t_arr - t0) <= float(uq_spec['event_half_window_s'])
    sigma = np.array([v.estimate_trace_snr(D[:, j], t_arr)['sigma'] for j in traces], dtype=float)
    qgrid = np.asarray(uq_spec['q_grid'], dtype=float)
    fgrid = np.asarray(uq_spec['frequency_factors'], dtype=float) * float(public['nominal_frequency'])
    huber_d = float(uq_spec['huber_delta'])
    w_prior = float(uq_spec['source_prior_weight'])
    s_prior = float(uq_spec['source_prior_log_scale'])
    nom_f = float(public['nominal_frequency'])

    c_unit = dict(public); c_unit['Q'] = [1.0, 1.0, 1.0, 1.0]
    travels = []
    for j in traces:
        a, b, _, _ = v.ray_absorption(times, float(offsets[j]), c_unit)
        travels.append(a + b)
    N = int(public['nfft']); dt = float(public['dt']); fref = float(public['fref'])
    f = np.arange(0, N // 2 + 1, dtype=float) / (N * dt)
    z_vec = -np.pi * f + 2j * f * np.log(np.maximum(f, np.finfo(float).tiny) / fref)
    z_vec[0] = 0.0
    phase_half = np.concatenate([[0.0], -1j * np.ones(N // 2 - 1), [0.0]])

    def build_templates_fast(jj, qq, ff):
        travel = travels[jj]
        B = np.zeros((nt, 6))
        for k in range(3):
            x = np.pi * ff * (t_arr - times[k])
            w = (1.0 - 2.0 * (x ** 2)) * np.exp(-(x ** 2))
            wf = np.fft.rfft(w, n=N)
            h = np.exp(z_vec * (travel[k] / qq))
            h[-1] = np.real(h[-1])
            u = wf * h
            B[:, 2 * k] = np.fft.irfft(u, n=N)[:nt]
            B[:, 2 * k + 1] = np.fft.irfft(u * phase_half, n=N)[:nt]
        return B

    def objective(qq, ff):
        loss = w_prior * ((math.log(ff / nom_f) / s_prior) ** 2)
        for jj, j in enumerate(traces):
            B = build_templates_fast(jj, qq, ff)
            sig = max(sigma[jj], np.finfo(float).tiny)
            X = B[mask, :] / sig
            y = D[mask, j] / sig
            weights = np.ones_like(y)
            coef = np.zeros(6)
            for _ in range(int(uq_spec['irls_iterations'])):
                WX = weights[:, None] * X
                G = X.T @ WX
                ridge = 1e-8 * max(np.trace(G) / 6.0, np.finfo(float).tiny)
                coef = sla.solve(G + ridge * np.eye(6), X.T @ (weights * y), assume_a='pos')
                r = y - X @ coef
                weights = np.minimum(1.0, huber_d / np.maximum(np.abs(r), np.finfo(float).tiny))
            r = np.abs(y - X @ coef)
            small = r <= huber_d
            loss += float(np.sum(0.5 * (r[small] ** 2)) + np.sum(huber_d * (r[~small] - 0.5 * huber_d)))
        return loss

    losses = np.zeros((len(qgrid), len(fgrid)))
    for iq, qq in enumerate(qgrid):
        for jf, ff in enumerate(fgrid):
            losses[iq, jf] = objective(qq, ff)
    idx = int(np.argmin(losses))
    iq, jf = np.unravel_index(idx, losses.shape)
    best = float(losses[iq, jf])
    q_best = float(qgrid[iq])
    f_best = float(fgrid[jf])

    def bounded(v_log):
        qq, ff = math.exp(v_log[0]), math.exp(v_log[1])
        if qq < qgrid[0] or qq > qgrid[-1] or ff < fgrid[0] or ff > fgrid[-1]:
            return 1e12 + float(np.sum(v_log ** 2))
        return objective(qq, ff)

    res = sopt.minimize(bounded, np.log([q_best, f_best]), method='Nelder-Mead',
                        options={'maxfev': 65, 'maxiter': 45, 'xatol': 1e-3, 'fatol': 1e-2})
    if np.isfinite(res.fun) and res.fun < best:
        q_best = float(math.exp(res.x[0]))
        f_best = float(math.exp(res.x[1]))
        best = float(res.fun)

    threshold = best + float(uq_spec['profile_loss_increment_per_sample']) * int(np.sum(mask)) * len(traces)
    qi, fi = np.where(losses <= threshold)
    qrange = [min(q_best, float(np.min(qgrid[qi]))), max(q_best, float(np.max(qgrid[qi])))] if len(qi) else [q_best, q_best]
    frange = [min(f_best, float(np.min(fgrid[fi]))), max(f_best, float(np.max(fgrid[fi])))] if len(fi) else [f_best, f_best]

    est = {
        'effective_Q': q_best,
        'estimated_frequency': f_best,
        'nominal_frequency': nom_f,
        'minimum_loss': best,
        'near_optimal_threshold': threshold,
        'q_sensitivity_range': qrange,
        'frequency_sensitivity_range': frange,
        'q_grid': qgrid.tolist(),
        'frequency_grid': fgrid.tolist(),
        'profile_losses': losses.tolist(),
        'estimated_sigma': sigma.tolist()
    }
    return enrich_calibration_uncertainty(est)

def solve_w_space_fast(sc, rt, z, u, V_l, perp, lam, gamma, HH, hh_diag_mean, r2_mat, extra_diag=None, rtol=1e-3, maxiter=25, w0_in=None, same_V=False, V_3d=None, VT_3d=None):
    n, nx = z.shape
    lam_s = sc * lam[None, :]
    mu_s = sc * gamma * float(np.median(lam))
    z2 = z ** 2
    d_mode_2d = z2 + lam_s * r2_mat
    d_mode = d_mode_2d.ravel(order='F')
    b_w = (z * u).ravel(order='F')
    if gamma == 0.0 or HH.nnz == 0:
        if extra_diag is not None and np.any(extra_diag != 0):
            W = np.zeros((n, nx))
            res_num = 0.0
            res_den = 0.0
            for j in range(nx):
                Vj = V_l[j]
                Nj = (Vj * d_mode_2d[:, j][None, :]) @ Vj.T + np.diag(sc * extra_diag[:, j])
                rhs_j = Vj @ (z[:, j] * u[:, j])
                sj = sla.solve(Nj, rhs_j, assume_a='pos')
                W[:, j] = Vj.T @ sj
                res_num += float(np.sum((Nj @ sj - rhs_j) ** 2))
                res_den += float(np.sum(rhs_j ** 2))
            normal = math.sqrt(res_num / max(res_den, 1e-30))
        else:
            W = (z * u) / d_mode_2d
            w_flat = W.ravel(order='F')
            normal = float(sla.norm(d_mode * w_flat - b_w) / max(sla.norm(b_w), 1e-30))
        rr = float(np.sum((z * W - u) ** 2) + np.sum(perp))
        return W, W.ravel(order='F'), {'normal_residual': normal, 'data_residual': float(math.sqrt(max(0.0, rr))), 'flag': 0, 'solver': 'SVD'}

    extra_mean = np.mean(extra_diag, axis=0) if extra_diag is not None else np.zeros(nx)
    inv_mode = 1.0 / (d_mode + ((mu_s * hh_diag_mean + sc * extra_mean)[None, :] * np.ones((n, 1))).ravel(order='F'))
    w0 = (b_w / d_mode) if w0_in is None else w0_in
    V0 = V_l[0]
    def matvec(w):
        w_flat = np.asarray(w, dtype=float).ravel()
        W = w_flat.reshape((n, nx), order='F')
        if same_V:
            S_m = V0 @ W
        elif V_3d is not None:
            S_m = (V_3d @ W.T[:, :, None]).squeeze(-1).T
        else:
            S_m = np.column_stack([V_l[j] @ W[:, j] for j in range(nx)])
        Hs = np.asarray(HH @ S_m.ravel(order='F')).ravel().reshape((n, nx), order='F')
        if extra_diag is not None:
            Hs = Hs + (sc / max(mu_s, 1e-30)) * (extra_diag * S_m)
        if same_V:
            VtHs = (V0.T @ Hs).ravel(order='F')
        elif VT_3d is not None:
            VtHs = (VT_3d @ Hs.T[:, :, None]).squeeze(-1).T.ravel(order='F')
        else:
            VtHs = np.column_stack([V_l[j].T @ Hs[:, j] for j in range(nx)]).ravel(order='F')
        return d_mode * w_flat + mu_s * VtHs

    def precond(r_vec):
        return inv_mode * np.asarray(r_vec, dtype=float).ravel()

    w, flag = spla.cg(
        spla.LinearOperator((n * nx, n * nx), matvec=matvec, dtype=float),
        b_w, x0=w0, rtol=rtol, atol=0.0, maxiter=maxiter,
        M=spla.LinearOperator((n * nx, n * nx), matvec=precond, dtype=float)
    )
    W = w.reshape((n, nx), order='F')
    rr = float(np.sum((z * W - u) ** 2) + np.sum(perp))
    normal = float(sla.norm(matvec(w) - b_w) / max(sla.norm(b_w), 1e-30))
    return W, w, {'normal_residual': normal, 'data_residual': float(math.sqrt(max(0.0, rr))), 'flag': int(flag), 'solver': 'PCG'}

def run_discrepancy_curve_and_solve(targets_dict, zD, zS, uD, uS, VD, VS, perpD, perpS, lam, gamma, HH, hh_diag_mean, r2_D, r2_S, extra_diag, v2_opts, same_V, VD_3d=None, VDT_3d=None, VS_3d=None, VST_3d=None):
    n, nx = zD.shape
    grid = list(map(float, v2_opts['scale_grid']))
    cache_r = {}
    def get_r(i):
        if i not in cache_r:
            cache_r[i] = solve_w_space_fast(grid[i], 'Direct', zD, uD, VD, perpD, lam, gamma, HH, hh_diag_mean, r2_D, extra_diag, rtol=2e-3, maxiter=18, same_V=same_V, V_3d=VD_3d, VT_3d=VDT_3d)
        return cache_r[i][2]['data_residual']

    results = {}
    for mname, target in targets_dict.items():
        r4 = get_r(4)
        if r4 <= target:
            c_i = 4
            while c_i + 1 < len(grid) and get_r(c_i + 1) <= target:
                c_i += 1
            below = c_i; above = c_i + 1 if c_i + 1 < len(grid) else None
        else:
            c_i = 4
            while c_i - 1 >= 0 and get_r(c_i - 1) >= target:
                c_i -= 1
            above = c_i; below = c_i - 1 if c_i - 1 >= 0 else None

        scales = [grid[i] for i in sorted(cache_r.keys())]
        residuals = [cache_r[i][2]['data_residual'] for i in sorted(cache_r.keys())]
        w_vecs = [cache_r[i][1] for i in sorted(cache_r.keys())]

        if below is not None and above is not None and below < above:
            lo = grid[below]; hi = grid[above]
            r_lo = cache_r[below][2]['data_residual']
            r_hi = cache_r[above][2]['data_residual']
            last_w = cache_r[below][1]
            for step in range(int(v2_opts['refinement_max_steps'])):
                if min(abs(r / target - 1.0) for r in residuals) <= float(v2_opts['residual_tolerance']):
                    break
                if step == 0 and r_hi > r_lo > 0 and target > 0:
                    frac = float(np.clip(math.log(target / r_lo) / max(math.log(r_hi / r_lo), 1e-12), 0.15, 0.85))
                    sc = lo * ((hi / lo) ** frac)
                else:
                    sc = math.sqrt(lo * hi)
                _, last_w, rec = solve_w_space_fast(sc, 'Direct', zD, uD, VD, perpD, lam, gamma, HH, hh_diag_mean, r2_D, extra_diag, rtol=2e-3, maxiter=18, w0_in=last_w, same_V=same_V, V_3d=VD_3d, VT_3d=VDT_3d)
                rv = rec['data_residual']
                scales.append(sc); residuals.append(rv); w_vecs.append(last_w)
                if rv < target:
                    lo = sc; r_lo = rv
                else:
                    hi = sc; r_hi = rv

        ix = int(np.argmin([abs(r - target) for r in residuals]))
        scale = scales[ix]
        W_d, w_d_vec, rec_d = solve_w_space_fast(scale, 'Direct', zD, uD, VD, perpD, lam, gamma, HH, hh_diag_mean, r2_D, extra_diag, rtol=1.1e-7, maxiter=200, w0_in=w_vecs[ix], same_V=same_V, V_3d=VD_3d, VT_3d=VDT_3d)
        W_s, _, rec_s = solve_w_space_fast(scale, 'Split', zS, uS, VS, perpS, lam, gamma, HH, hh_diag_mean, r2_S, extra_diag, rtol=1.1e-7, maxiter=200, w0_in=w_d_vec, same_V=same_V, V_3d=VS_3d, VT_3d=VST_3d)
        if same_V:
            S_d = np.ascontiguousarray(VD[0] @ W_d)
            S_s = np.ascontiguousarray(VS[0] @ W_s)
        else:
            S_d = np.ascontiguousarray((VD_3d @ W_d.T[:, :, None]).squeeze(-1).T)
            S_s = np.ascontiguousarray((VS_3d @ W_s.T[:, :, None]).squeeze(-1).T)
        endpoint = int(scale == grid[0] or scale == grid[-1])
        mismatch = float(abs(residuals[ix] / target - 1.0))
        results[mname] = {
            'direct': S_d, 'split': S_s, 'scale': scale, 'target': target,
            'endpoint': endpoint, 'residual_mismatch': mismatch,
            'direct_record': rec_d, 'split_record': rec_s
        }
    return results

def evaluate_all_branches_on_gather(D, models, domain, cm, R_avo, o_struct, v2_opts, p1_hyp):
    n, nx = D.shape
    same_V = (domain == 'Section')
    sigma = np.array([v.estimate_trace_snr(D[:, j], cm['t'])['sigma'] for j in range(nx)], dtype=float)
    lam = np.array([v.discrepancy_apply(models[j]['direct'], D[:, j], cm['lambda_grid'], math.sqrt(n) * sigma[j])[1] for j in range(nx)], dtype=float)
    pilot = np.ascontiguousarray(np.column_stack([
        models[j]['direct']['V'] @ ((models[j]['direct']['z'] / (models[j]['direct']['z'] ** 2 + 4.0 * lam[j])) * (models[j]['direct']['U'].T @ D[:, j]))
        for j in range(nx)
    ]))

    zD = np.column_stack([models[j]['direct']['z'] for j in range(nx)])
    zS = np.column_stack([models[j]['split']['z'] for j in range(nx)])
    VD = [models[j]['direct']['V'] for j in range(nx)]
    VS = [models[j]['split']['V'] for j in range(nx)]
    VD_3d = None if same_V else np.ascontiguousarray(np.stack(VD, axis=0))
    VDT_3d = None if same_V else np.ascontiguousarray(np.transpose(VD_3d, (0, 2, 1)))
    VS_3d = None if same_V else np.ascontiguousarray(np.stack(VS, axis=0))
    VST_3d = None if same_V else np.ascontiguousarray(np.transpose(VS_3d, (0, 2, 1)))
    ones_r2 = np.ones((n, nx))
    base_norm = math.sqrt(n * float(np.sum(sigma ** 2)))

    Yw_orig = np.ascontiguousarray(np.column_stack([models[j]['K'] @ D[:, j] for j in range(nx)]))
    uD_orig = np.ascontiguousarray(np.column_stack([models[j]['direct']['U'].T @ D[:, j] for j in range(nx)]))
    uS_orig = np.ascontiguousarray(np.column_stack([models[j]['split']['U'].T @ Yw_orig[:, j] for j in range(nx)]))
    perpD_orig = np.maximum(0.0, np.sum(D ** 2, axis=0) - np.sum(uD_orig ** 2, axis=0))
    perpS_orig = np.maximum(0.0, np.sum(Yw_orig ** 2, axis=0) - np.sum(uS_orig ** 2, axis=0))

    outputs = {}
    H_zero = sp.csc_matrix((n * nx, n * nx))

    res_base = run_discrepancy_curve_and_solve(
        {'Baseline': 1.0 * base_norm},
        zD, zS, uD_orig, uS_orig, VD, VS, perpD_orig, perpS_orig,
        lam, 0.0, H_zero, np.zeros(nx), ones_r2, ones_r2, None, v2_opts, same_V, VD_3d, VDT_3d, VS_3d, VST_3d
    )['Baseline']
    outputs['Baseline'] = (res_base, {'H': H_zero, 'sigma': sigma, 'lambda': lam, 'pilot': pilot},
                           {'name': 'Baseline', 'gamma': 0.0, 'pilot_multiplier': 4.0, 'noise_target_multiplier': 1.0})

    H0_orig = v.learn_structure_operator(pilot, cm, o_struct)
    H_froz = (H0_orig @ R_avo).tocsc() if domain == 'Offset' else H0_orig.tocsc()
    HH_froz = (H_froz.T @ H_froz).tocsc()
    hh_diag_froz = HH_froz.diagonal()
    hh_mean_froz = np.array([float(np.mean(hh_diag_froz[j*n:(j+1)*n])) for j in range(nx)])
    res_fs = run_discrepancy_curve_and_solve(
        {'Frozen': 1.0 * base_norm, 'StrongerDamping': 1.15 * base_norm},
        zD, zS, uD_orig, uS_orig, VD, VS, perpD_orig, perpS_orig,
        lam, 10.0, HH_froz, hh_mean_froz, ones_r2, ones_r2, None, v2_opts, same_V, VD_3d, VDT_3d, VS_3d, VST_3d
    )
    outputs['Frozen'] = (res_fs['Frozen'], {'H': H_froz, 'sigma': sigma, 'lambda': lam, 'pilot': pilot},
                         {'name': 'Frozen', 'gamma': 10.0, 'pilot_multiplier': 4.0, 'noise_target_multiplier': 1.0})
    outputs['StrongerDamping'] = (res_fs['StrongerDamping'], {'H': H_froz, 'sigma': sigma, 'lambda': lam, 'pilot': pilot},
                                  {'name': 'StrongerDamping', 'gamma': 10.0, 'pilot_multiplier': 4.0, 'noise_target_multiplier': 1.15})

    b_coh = np.zeros_like(D)
    extra_diag = None
    tail_mask = np.zeros((n, nx), dtype=bool)
    if domain == 'Offset':
        pilot_ortho = (R_avo @ pilot.ravel(order='F')).reshape((n, nx), order='F')
        D_ortho = np.column_stack([models[j]['G'] @ pilot_ortho[:, j] for j in range(nx)])
        expected_ortho = (nx - 2) / nx * n * float(np.sum(sigma ** 2))
        ortho_ratio = float(np.sum(D_ortho ** 2) / max(expected_ortho, 1e-30))
        if ortho_ratio >= float(p1_hyp['coherent_ortho_ratio_gate']):
            t_idx = np.arange(n); lags = p1_hyp['coherent_prediction_lags']
            for j in range(nx):
                neighbors = [jj for jj in (j - 1, j + 1) if 0 <= jj < nx]
                B_mat = np.column_stack([D_ortho[np.clip(t_idx + lag, 0, n - 1), jj] for jj in neighbors for lag in lags])
                p_col = B_mat.shape[1]; ridge = float(p1_hyp['coherent_prediction_ridge']) * max(np.sum(B_mat ** 2) / p_col, 1e-30)
                coef = sla.solve(B_mat.T @ B_mat + ridge * np.eye(p_col), B_mat.T @ D_ortho[:, j], assume_a='pos')
                pred = B_mat @ coef
                r2_c = 1.0 - np.sum((D_ortho[:, j] - pred) ** 2) / max(np.sum(D_ortho[:, j] ** 2), 1e-30)
                b_coh[:, j] = np.clip((r2_c - float(p1_hyp['coherent_r2_gate'])) / float(p1_hyp['coherent_r2_span']), 0.0, 1.0) * pred
    else:
        loc_pow = smooth_box(pilot ** 2, 19)
        lat_med = np.median(loc_pow, axis=1, keepdims=True)
        sigma_p2 = np.array([(sigma[j] ** 2) * np.mean((models[j]['direct']['z'] / (models[j]['direct']['z'] ** 2 + 4.0 * lam[j])) ** 2) for j in range(nx)])[None, :]
        tail_mask = (loc_pow < float(p1_hyp['section_tail_power_factor']) * sigma_p2) & (loc_pow < float(p1_hyp['section_tail_lateral_ratio']) * lat_med)
        extra_diag = np.zeros((n, nx))
        for j in range(nx):
            extra_diag[:, j] = np.where(tail_mask[:, j], float(p1_hyp['section_tail_diag_lambda_mult']) * lam[j], 0.0)

    D_work = D - b_coh
    if np.any(b_coh != 0):
        pilot_w = np.ascontiguousarray(np.column_stack([
            models[j]['direct']['V'] @ ((models[j]['direct']['z'] / (models[j]['direct']['z'] ** 2 + 4.0 * lam[j])) * (models[j]['direct']['U'].T @ D_work[:, j]))
            for j in range(nx)
        ]))
        Yw_w = np.ascontiguousarray(np.column_stack([models[j]['K'] @ D_work[:, j] for j in range(nx)]))
        uD_w = np.ascontiguousarray(np.column_stack([models[j]['direct']['U'].T @ D_work[:, j] for j in range(nx)]))
        uS_w = np.ascontiguousarray(np.column_stack([models[j]['split']['U'].T @ Yw_w[:, j] for j in range(nx)]))
        perpD_w = np.maximum(0.0, np.sum(D_work ** 2, axis=0) - np.sum(uD_w ** 2, axis=0))
        perpS_w = np.maximum(0.0, np.sum(Yw_w ** 2, axis=0) - np.sum(uS_w ** 2, axis=0))
        H0_w = v.learn_structure_operator(pilot_w, cm, o_struct)
    else:
        pilot_w = pilot; Yw_w = Yw_orig; uD_w = uD_orig; uS_w = uS_orig; perpD_w = perpD_orig; perpS_w = perpS_orig; H0_w = H0_orig

    # Directly incorporate calibrated (Q, fpeak) posterior uncertainty (cm['q_uncertainty'], cm['f_uncertainty'], cm['log_qf_corr'])
    fpeak = max(float(cm['fpeak']), 10.0); q_eff = float(cm['Q'][0])
    low_q_weight = math.exp(-max(0.0, q_eff - float(p1_hyp['low_q_reference'])) / float(p1_hyp['low_q_decay_scale']))
    u_q = float(cm.get('q_uncertainty', 0.12))
    u_f_cal = float(cm.get('f_uncertainty', 0.12))
    rho_qf = float(cm.get('log_qf_corr', -0.75))
    uq_norm = float(np.clip(u_q / 0.12, 0.75, 1.35))
    uf_norm = float(np.clip(u_f_cal / 0.14, 0.75, 1.35))
    rho_abs = float(np.clip(abs(rho_qf), 0.0, 0.95))

    r2_noise_D = np.zeros((n, nx)); r2_noise_S = np.zeros((n, nx))
    r2_q_D = np.zeros((n, nx)); r2_q_S = np.zeros((n, nx))
    for j in range(nx):
        for rk, arr_n, arr_q in (('direct', r2_noise_D, r2_q_D), ('split', r2_noise_S, r2_q_S)):
            F = models[j][rk]; z = F['z']; f_m = F['f_mode']
            u_f = f_m / fpeak
            s_prior = (u_f ** 2) * np.exp(1.0 - np.minimum(u_f ** 2, 20.0))
            hf_n = 0.6 * (np.maximum(0.0, u_f - float(p1_hyp['hf_cutoff_factor_noise'])) ** 2)
            arr_n[:, j] = (1.0 + hf_n) / (1.0 + 0.06 * s_prior)
            j_q = np.log(1.0 / np.maximum(z, 1e-4))
            j_f = np.maximum(0.0, u_f - float(p1_hyp['hf_cutoff_factor_q']))
            q_sens = (float(p1_hyp['q_sensitivity_scale']) * low_q_weight * uq_norm * j_q) ** 2
            hf_q = (2.4 * low_q_weight * uq_norm + 0.38 * uf_norm) * (j_f ** 2) + 0.05 * rho_abs * low_q_weight * uq_norm * uf_norm * (j_q * j_f)
            pass_boost = (float(p1_hyp['passband_boost_lowq']) * low_q_weight + float(p1_hyp['passband_boost_base'])) / (uf_norm ** 1.2)
            arr_q[:, j] = (1.0 + hf_q + q_sens) / (1.0 + pass_boost * s_prior)

    uq_eta = 0.75 + 0.25 * uq_norm
    eta_noise = float(p1_hyp['eta_noise_only_offset']) if domain == 'Offset' else float(p1_hyp['eta_noise_only_section'])
    eta_q = (float(p1_hyp['eta_offset_base']) + float(p1_hyp['eta_offset_lowq_add']) * low_q_weight * uq_eta) if domain == 'Offset' else (float(p1_hyp['eta_section_base']) + float(p1_hyp['eta_section_lowq_add']) * low_q_weight * uq_eta + 0.025 * max(0.0, uf_norm - 1.0))

    if domain == 'Offset':
        H_comb = sp.vstack([H0_w @ R_avo, float(p1_hyp['offset_avo_ortho_weight_noise']) * R_avo]).tocsc()
        H_qonly = sp.vstack([H0_orig @ R_avo, float(p1_hyp['offset_avo_ortho_weight_qonly']) * R_avo]).tocsc()
    else:
        if np.any(tail_mask):
            row_scale = np.where(tail_mask, float(p1_hyp['section_tail_h_scale']), 1.0).ravel(order='F')
            H_comb = H0_w.multiply(row_scale[:, None]).tocsc()
            H_qonly = H0_orig.multiply(row_scale[:, None]).tocsc()
        else:
            H_comb = H0_w.tocsc(); H_qonly = H0_orig.tocsc()

    HH_comb = (H_comb.T @ H_comb).tocsc()
    hh_diag_comb = HH_comb.diagonal()
    hh_mean_comb = np.array([float(np.mean(hh_diag_comb[j*n:(j+1)*n])) for j in range(nx)])

    res_no = run_discrepancy_curve_and_solve(
        {'NoiseUncertaintyOnly': eta_noise * base_norm},
        zD, zS, uD_w, uS_w, VD, VS, perpD_w, perpS_w,
        lam, 10.0, HH_comb, hh_mean_comb, r2_noise_D, r2_noise_S, extra_diag, v2_opts, same_V, VD_3d, VDT_3d, VS_3d, VST_3d
    )['NoiseUncertaintyOnly']
    outputs['NoiseUncertaintyOnly'] = (res_no, {'H': H_comb, 'sigma': sigma, 'lambda': lam, 'pilot': pilot_w},
                                       {'name': 'NoiseUncertaintyOnly', 'gamma': 10.0, 'pilot_multiplier': 4.0, 'noise_target_multiplier': eta_noise})

    res_cu = run_discrepancy_curve_and_solve(
        {'CombinedUncertainty': eta_q * base_norm},
        zD, zS, uD_w, uS_w, VD, VS, perpD_w, perpS_w,
        lam, 10.0, HH_comb, hh_mean_comb, r2_q_D, r2_q_S, extra_diag, v2_opts, same_V, VD_3d, VDT_3d, VS_3d, VST_3d
    )['CombinedUncertainty']
    outputs['CombinedUncertainty'] = (res_cu, {'H': H_comb, 'sigma': sigma, 'lambda': lam, 'pilot': pilot_w},
                                      {'name': 'CombinedUncertainty', 'gamma': 10.0, 'pilot_multiplier': 4.0, 'noise_target_multiplier': eta_q})

    res_cuns = run_discrepancy_curve_and_solve(
        {'CombinedUncertainty_NoStructure': eta_q * base_norm},
        zD, zS, uD_w, uS_w, VD, VS, perpD_w, perpS_w,
        lam, 0.0, H_zero, np.zeros(nx), r2_q_D, r2_q_S, None, v2_opts, same_V, VD_3d, VDT_3d, VS_3d, VST_3d
    )['CombinedUncertainty_NoStructure']
    outputs['CombinedUncertainty_NoStructure'] = (res_cuns, {'H': H_zero, 'sigma': sigma, 'lambda': lam, 'pilot': pilot_w},
                                                  {'name': 'CombinedUncertainty_NoStructure', 'gamma': 0.0, 'pilot_multiplier': 4.0, 'noise_target_multiplier': eta_q})

    HH_qonly = (H_qonly.T @ H_qonly).tocsc()
    hh_diag_qonly = HH_qonly.diagonal()
    hh_mean_qonly = np.array([float(np.mean(hh_diag_qonly[j*n:(j+1)*n])) for j in range(nx)])
    res_qo = run_discrepancy_curve_and_solve(
        {'QModelErrorOnly': eta_q * base_norm},
        zD, zS, uD_orig, uS_orig, VD, VS, perpD_orig, perpS_orig,
        lam, 10.0, HH_qonly, hh_mean_qonly, r2_q_D, r2_q_S, extra_diag, v2_opts, same_V, VD_3d, VDT_3d, VS_3d, VST_3d
    )['QModelErrorOnly']
    outputs['QModelErrorOnly'] = (res_qo, {'H': H_qonly, 'sigma': sigma, 'lambda': lam, 'pilot': pilot},
                                  {'name': 'QModelErrorOnly', 'gamma': 10.0, 'pilot_multiplier': 4.0, 'noise_target_multiplier': eta_q})

    return outputs
