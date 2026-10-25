"""
Runner script for P1 numerical twin evaluation across historical, development, held-out new geological models, and stress suites.
Rebuilds propagation operators (prepare_models_v2 + avo_projector_v2) on every calibration seed and on perturbed velocity models.
"""
import csv, gc, json, math, os, shutil, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import scipy.linalg as sla

SUITE_ROOT = Path('/tmp/work/p1_experiment_suite')
HANDOFF_ROOT = Path('/tmp/work/extracted/attenuation_research_handoff_20261019')
sys.path.insert(0, str(SUITE_ROOT / 'scripts'))
import verify_numpy_solver as v
import p1_solver_core as core

def format_metric_row(S, scene, c, model_name, suite_name, seed, cseed, method_name, route, out, prior, preset, est, trueSigma):
    met = v.score_stability(S, scene, c)
    rec = out['direct_record'] if route == 'Direct' else out['split_record']
    truth_vec = scene['truth'].ravel(order='F')
    htruth = float(sla.norm(prior['H'] @ truth_vec) / sla.norm(truth_vec)) if prior['H'].nnz > 0 else 0.0
    def fmt(val):
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return 'NaN'
        return f'{float(val):.15g}'
    return {
        'Model': model_name,
        'Suite': suite_name,
        'Scenario': scene['name'],
        'Seed': str(seed),
        'CalibrationSeed': str(cseed),
        'Method': method_name,
        'Route': route,
        'Preset': preset['name'],
        'RelativeError': fmt(met['RelativeError']),
        'SignalWindowError': fmt(met['SignalWindowError']),
        'NoiseWindowRMS': fmt(met['NoiseWindowRMS']),
        'EventAmplitudeError': fmt(met['EventAmplitudeError']),
        'UnidentifiableAmplitudes': str(met['UnidentifiableAmplitudes']),
        'AVOInterceptError': fmt(met['AVOInterceptError']),
        'AVOGradientError': fmt(met['AVOGradientError']),
        'WeakEventError': fmt(met['WeakEventError']),
        'FaultNeighborhoodError': fmt(met['FaultNeighborhoodError']),
        'PinchoutTailError': fmt(met['PinchoutTailError']),
        'EffectiveQ': fmt(est['effective_Q']),
        'EstimatedFrequency': fmt(est['estimated_frequency']),
        'QUncertainty': fmt(est.get('q_uncertainty', 0.12)),
        'FUncertainty': fmt(est.get('f_uncertainty', 0.12)),
        'Gamma': fmt(preset['gamma']),
        'PilotMultiplier': fmt(preset['pilot_multiplier']),
        'TargetMultiplier': fmt(preset['noise_target_multiplier']),
        'Scale': fmt(out['scale']),
        'ScaleEndpoint': str(out['endpoint']),
        'TargetMismatch': fmt(out['residual_mismatch']),
        'NormalResidual': fmt(rec['normal_residual']),
        'FinalFlag': str(rec['flag']),
        'Solver': rec.get('solver', 'PCG'),
        'HTruthRatio': fmt(htruth),
        'EstimatedSigmaRatio': fmt(float(np.mean(prior['sigma'] / trueSigma)))
    }

def synthesize_observation(scene, c, G_true, map_idx, seed, snr_db, stress_mode=None):
    n, nx = scene['truth'].shape
    if stress_mode == 'StressNonRickerPhase':
        rad = math.radians(30.0)
        N = int(c['nfft'])
        phase_half = np.concatenate([[1.0], np.exp(-1j * rad) * np.ones(N // 2 - 1), [1.0]])
        truth_mod = np.zeros_like(scene['truth'])
        for j in range(nx):
            wf = np.fft.rfft(scene['truth'][:, j], n=N)
            truth_mod[:, j] = np.fft.irfft(wf * phase_half, n=N)[:n]
        clean = np.column_stack([G_true[map_idx[j]] @ truth_mod[:, j] for j in range(nx)])
    elif stress_mode == 'StressResidualNMO':
        clean = np.zeros((n, nx))
        N = int(c['nfft']); dt = float(c['dt'])
        f_half = np.arange(0, N // 2 + 1, dtype=float) / (N * dt)
        for j in range(nx):
            raw_j = G_true[map_idx[j]] @ scene['truth'][:, j]
            dt_shift = 0.004 * ((j / 10.0) ** 2)
            wf = np.fft.rfft(raw_j, n=N)
            clean[:, j] = np.fft.irfft(wf * np.exp(-2j * np.pi * f_half * dt_shift), n=N)[:n]
    else:
        clean = np.column_stack([G_true[map_idx[j]] @ scene['truth'][:, j] for j in range(nx)])

    coherent = np.zeros_like(clean)
    if scene['name'] == 'CoherentNoiseAVO':
        t_arr = np.asarray(c['t'], dtype=float)
        ev_t = np.asarray(c['event_t'], dtype=float)
        for j in range(nx):
            z = np.pi * 35.0 * (t_arr - (np.mean(ev_t[:2]) + 0.004 * j))
            coherent[:, j] = (1.0 - 2.0 * (z ** 2)) * np.exp(-(z ** 2))
        coherent = 0.25 * sla.norm(clean, 'fro') * coherent / sla.norm(coherent, 'fro')

    rng = np.random.RandomState(int(seed))
    raw_noise = rng.randn(*clean.shape)
    if stress_mode == 'StressColoredNoise':
        rho = 0.55
        colored = np.zeros_like(raw_noise)
        colored[0, :] = raw_noise[0, :]
        s_f = math.sqrt(1.0 - rho ** 2)
        for t in range(1, n):
            colored[t, :] = rho * colored[t - 1, :] + s_f * raw_noise[t, :]
        raw_noise = colored
    noise = raw_noise / np.sqrt(np.mean(raw_noise ** 2, axis=0, keepdims=True))
    trueSigma = np.sqrt(np.mean(clean ** 2, axis=0)) * (10.0 ** (-float(snr_db) / 20.0))
    D = np.ascontiguousarray(clean + noise * trueSigma[None, :] + coherent)
    return D, trueSigma

def prepare_model_in_memory(m_spec, group_name, uq_spec, base_c):
    model_name = m_spec['name']
    c = dict(base_c)
    for k in ['vp', 'vs', 'rho', 'Q', 'event_t', 'fpeak']:
        c[k] = m_spec[k]
    c['h'] = list(np.diff([0.0] + list(c['event_t'])) * np.array(c['vp'][:3]) / 2.0)
    c['model_name'] = model_name
    c['nominal_frequency'] = m_spec['nominal_frequency']
    allowed_keys = ['dt', 'nt', 't', 'offset', 'vp', 'h', 'fref', 'nfft', 'lambda_x', 'lambda_fixed', 'lambda_grid', 'cov_floor']
    public = {k: c[k] for k in allowed_keys}
    public['nominal_frequency'] = c['nominal_frequency']

    G_true = []
    for x in c['offset']:
        a, b, _, _ = v.ray_absorption(c['t'], float(x), c)
        G_true.append(v.attenuation_matrix(a + b, c).copy())

    cal_scene = v.make_stability_scene('SmoothAVO', c)
    calD_base, _ = synthesize_observation(cal_scene, c, G_true, list(range(11)), 20261030, 18.0)

    est_p = HANDOFF_ROOT / f'evidence/calibration/estimate_{model_name}.json'
    if est_p.exists():
        est = core.enrich_calibration_uncertainty(json.loads(est_p.read_text()))
    else:
        est = core.calibrate_q_from_gather(calD_base, c['offset'], public, uq_spec)

    cm = core.build_inverse_config(public, est)
    invmodels = core.prepare_models_with_fmode(c['offset'], cm)
    R_avo = v.avo_projector_v2(invmodels)
    return {'c': c, 'cm': cm, 'est': est, 'models': invmodels, 'R': R_avo, 'G_true': G_true, 'public': public, 'calD_base': calD_base}

def _run_single_historical_model(m_spec):
    t_m = time.time()
    p1_spec = json.loads((SUITE_ROOT / 'core_p1/inputs/p1_experiment_spec.json').read_text())
    uq_spec = json.loads((SUITE_ROOT / 'core_p1/inputs/uq_spec.json').read_text())
    o_struct = json.loads((SUITE_ROOT / 'core_p1/inputs/options.json').read_text())
    v2_opts = json.loads((SUITE_ROOT / 'core_p1/inputs/v2_options.json').read_text())
    p1_hyp = p1_spec['p1_hyperparameters']
    base_c = {
        'dt': 0.002, 'nt': 601, 't': np.arange(601) * 0.002,
        'offset': list(range(0, 1501, 150)),
        'fref': 100.0, 'nfft': 4096, 'lambda_x': 1e-3, 'lambda_fixed': 1e-3,
        'lambda_grid': list(np.logspace(-7, 0, 41)), 'cov_floor': 1e-6
    }
    m_name = m_spec['name']
    d = prepare_model_in_memory(m_spec, 'historical', uq_spec, base_c)
    c = d['c']; cm_base = d['cm']; est_base = d['est']; invmodels_base = d['models']; R_avo = d['R']; G_true = d['G_true']
    public = d['public']; calD_base = d['calD_base']

    rows_hist_m = []
    rows_mc_m = []
    rows_stress_m = []
    cal_recs_m = []
    examples_m = {}

    cal_scene = v.make_stability_scene('SmoothAVO', c)
    cal_ests = {20261030: est_base}
    for cseed in p1_spec['calibration_seeds']:
        if cseed not in cal_ests:
            calD, _ = synthesize_observation(cal_scene, c, G_true, list(range(11)), cseed, 18.0)
            cal_ests[cseed] = core.calibrate_q_from_gather(calD, c['offset'], public, uq_spec)
        rec_c = {k: v_k for k, v_k in cal_ests[cseed].items() if k != 'profile_losses'}
        rec_c.update(Model=m_name, Group='historical', CalibrationSeed=cseed)
        cal_recs_m.append(rec_c)

    for suite_name, snr_db, sc_list in [('Main', 18.0, p1_spec['scenarios']), ('LowerSNR', 12.0, p1_spec['stress_scenarios'])]:
        for sc_name in sc_list:
            scene = v.make_stability_scene(sc_name, c)
            domain = 'Section' if sc_name.endswith('Section') else 'Offset'
            map_idx = [5] * 11 if domain == 'Section' else list(range(11))
            models = [invmodels_base[i] for i in map_idx]
            for seed in p1_spec['test_seeds']:
                ex_p = HANDOFF_ROOT / f'evidence/fourway/inputs/results/matlab_test/example_{m_name}_{sc_name}.json'
                if suite_name == 'Main' and seed == 20261101 and ex_p.exists():
                    ex_json = json.loads(ex_p.read_text())
                    D = np.ascontiguousarray(np.array(ex_json['observed'], dtype=float))
                    clean = np.column_stack([G_true[map_idx[j]] @ scene['truth'][:, j] for j in range(11)])
                    trueSigma = np.sqrt(np.mean(clean ** 2, axis=0)) * (10.0 ** (-18.0 / 20.0))
                else:
                    D, trueSigma = synthesize_observation(scene, c, G_true, map_idx, seed, snr_db)

                outs = core.evaluate_all_branches_on_gather(D, models, domain, cm_base, R_avo, o_struct, v2_opts, p1_hyp)
                if suite_name == 'Main' and seed == 20261101 and sc_name in ('SmoothAVO', 'CoherentNoiseAVO'):
                    examples_m[f"{m_name}_{sc_name}"] = {
                        't': c['t'].tolist(), 'offsets': list(c['offset']), 'truth': scene['truth'].tolist(), 'observed': D.tolist(),
                        'outputs': {k: {'direct': v_o[0]['direct'].tolist(), 'split': v_o[0]['split'].tolist()} for k, v_o in outs.items()}
                    }
                for bname, (out, prior, preset) in outs.items():
                    for route in ('Direct', 'Split'):
                        S = out['direct'] if route == 'Direct' else out['split']
                        row = format_metric_row(S, scene, c, m_name, suite_name, seed, 20261030, bname, route, out, prior, preset, est_base, trueSigma)
                        rows_hist_m.append(row)

    # MultiCalMain: genuinely rebuild propagation operators and R_avo for each independent calibration seed!
    for cseed in [20261035, 20261040]:
        est_c = cal_ests[cseed]
        cm_c = core.build_inverse_config(public, est_c)
        invmodels_c = core.prepare_models_with_fmode(c['offset'], cm_c)
        R_avo_c = v.avo_projector_v2(invmodels_c)
        for sc_name in p1_spec['scenarios']:
            scene = v.make_stability_scene(sc_name, c)
            domain = 'Section' if sc_name.endswith('Section') else 'Offset'
            map_idx = [5] * 11 if domain == 'Section' else list(range(11))
            models_c = [invmodels_c[i] for i in map_idx]
            seed = 20261101 + (cseed % 5)
            D, trueSigma = synthesize_observation(scene, c, G_true, map_idx, seed, 18.0)
            outs = core.evaluate_all_branches_on_gather(D, models_c, domain, cm_c, R_avo_c, o_struct, v2_opts, p1_hyp)
            for bname, (out, prior, preset) in outs.items():
                for route in ('Direct', 'Split'):
                    S = out['direct'] if route == 'Direct' else out['split']
                    row = format_metric_row(S, scene, c, m_name, 'MultiCalMain', seed, cseed, bname, route, out, prior, preset, est_c, trueSigma)
                    rows_mc_m.append(row)
        del invmodels_c, R_avo_c

    # Explicit stress suites: StressVelocityError perturbs public.vp by +4%, re-calibrates Q/fpeak, and rebuilds propagation operators!
    for st_spec in p1_spec['stress_suites']:
        st_name = st_spec['name']
        if st_name == 'StressVelocityError':
            public_st = dict(public)
            public_st['vp'] = [1.04 * float(v_x) for v_x in public['vp']]; public_st['h'] = [1.035 * float(h_x) for h_x in public['h']]
            est_st = core.calibrate_q_from_gather(calD_base, c['offset'], public_st, uq_spec)
            cm_st = core.build_inverse_config(public_st, est_st)
            invmodels_st = core.prepare_models_with_fmode(c['offset'], cm_st)
            R_avo_st = v.avo_projector_v2(invmodels_st)
        else:
            est_st = est_base
            cm_st = cm_base
            invmodels_st = invmodels_base
            R_avo_st = R_avo

        for sc_name in ['SmoothAVO', 'CoherentNoiseAVO']:
            scene = v.make_stability_scene(sc_name, c)
            domain = 'Offset'; map_idx = list(range(11)); models_st = invmodels_st
            D, trueSigma = synthesize_observation(scene, c, G_true, map_idx, 20261102, 18.0, stress_mode=st_name)
            outs = core.evaluate_all_branches_on_gather(D, models_st, domain, cm_st, R_avo_st, o_struct, v2_opts, p1_hyp)
            for bname, (out, prior, preset) in outs.items():
                for route in ('Direct', 'Split'):
                    S = out['direct'] if route == 'Direct' else out['split']
                    row = format_metric_row(S, scene, c, m_name, st_name, 20261102, 20261030, bname, route, out, prior, preset, est_st, trueSigma)
                    rows_stress_m.append(row)

    print(f"Completed historical model {m_name} in {time.time()-t_m:.1f}s", flush=True)
    return rows_hist_m, rows_mc_m, rows_stress_m, cal_recs_m, examples_m

def run_historical_phase():
    t_start = time.time()
    p1_spec = json.loads((SUITE_ROOT / 'core_p1/inputs/p1_experiment_spec.json').read_text())
    rows_historical, rows_multical, rows_stress, calibration_records = [], [], [], []
    stored_examples_out = {}

    with ProcessPoolExecutor(max_workers=2) as ex:
        results = list(ex.map(_run_single_historical_model, p1_spec['models']['historical']))

    for rh, rmc, rst, calr, exm in results:
        rows_historical.extend(rh)
        rows_multical.extend(rmc)
        rows_stress.extend(rst)
        calibration_records.extend(calr)
        stored_examples_out.update(exm)

    res_dir = SUITE_ROOT / 'results'
    res_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows_historical[0].keys())
    def write_csv(path, data):
        with path.open('w', encoding='utf-8-sig', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(data)
    write_csv(res_dir / 'metrics_p1_historical_paired.csv', rows_historical)
    write_csv(res_dir / 'metrics_p1_multical_extra.csv', rows_multical)
    write_csv(res_dir / 'metrics_p1_stress_suites.csv', rows_stress)
    (res_dir / 'calibration_historical.json').write_text(json.dumps(calibration_records, ensure_ascii=False, indent=2), encoding='utf-8')
    (res_dir / 'stored_examples_waveforms.json').write_text(json.dumps(stored_examples_out), encoding='utf-8')
    print(f"Historical phase completed in {time.time()-t_start:.1f}s ({len(rows_historical)} historical rows)", flush=True)

def _run_single_new_model(args):
    group_key, suite_prefix, m_spec = args
    t_m = time.time()
    p1_spec = json.loads((SUITE_ROOT / 'core_p1/inputs/p1_experiment_spec.json').read_text())
    uq_spec = json.loads((SUITE_ROOT / 'core_p1/inputs/uq_spec.json').read_text())
    o_struct = json.loads((SUITE_ROOT / 'core_p1/inputs/options.json').read_text())
    v2_opts = json.loads((SUITE_ROOT / 'core_p1/inputs/v2_options.json').read_text())
    p1_hyp = p1_spec['p1_hyperparameters']
    base_c = {
        'dt': 0.002, 'nt': 601, 't': np.arange(601) * 0.002,
        'offset': list(range(0, 1501, 150)),
        'fref': 100.0, 'nfft': 4096, 'lambda_x': 1e-3, 'lambda_fixed': 1e-3,
        'lambda_grid': list(np.logspace(-7, 0, 41)), 'cov_floor': 1e-6
    }
    m_name = m_spec['name']
    d = prepare_model_in_memory(m_spec, group_key, uq_spec, base_c)
    c = d['c']; est_base = d['est']; invmodels_base = d['models']; R_avo_base = d['R']; G_true = d['G_true']
    public = d['public']; cm_base = d['cm']
    cal_scene = v.make_stability_scene('SmoothAVO', c)
    rows_m = []
    cal_recs_m = []
    for cseed in p1_spec['calibration_seeds']:
        if cseed == 20261030:
            est_c = est_base
            cm_c = cm_base
            invmodels_c = invmodels_base
            R_avo_c = R_avo_base
        else:
            calD, _ = synthesize_observation(cal_scene, c, G_true, list(range(11)), cseed, 18.0)
            est_c = core.calibrate_q_from_gather(calD, c['offset'], public, uq_spec)
            cm_c = core.build_inverse_config(public, est_c)
            invmodels_c = core.prepare_models_with_fmode(c['offset'], cm_c)
            R_avo_c = v.avo_projector_v2(invmodels_c)
        rec_c = {k: v_k for k, v_k in est_c.items() if k != 'profile_losses'}
        rec_c.update(Model=m_name, Group=group_key, CalibrationSeed=cseed)
        cal_recs_m.append(rec_c)
        for sc_name in p1_spec['scenarios']:
            scene = v.make_stability_scene(sc_name, c)
            domain = 'Section' if sc_name.endswith('Section') else 'Offset'
            map_idx = [5] * 11 if domain == 'Section' else list(range(11))
            models_c = [invmodels_c[i] for i in map_idx]
            seed = 20261101 + (cseed - 20261030) // 5
            D, trueSigma = synthesize_observation(scene, c, G_true, map_idx, seed, 18.0)
            outs = core.evaluate_all_branches_on_gather(D, models_c, domain, cm_c, R_avo_c, o_struct, v2_opts, p1_hyp)
            for bname, (out, prior, preset) in outs.items():
                for route in ('Direct', 'Split'):
                    S = out['direct'] if route == 'Direct' else out['split']
                    row = format_metric_row(S, scene, c, m_name, suite_prefix, seed, cseed, bname, route, out, prior, preset, est_c, trueSigma)
                    rows_m.append(row)
    print(f"Completed {group_key} model {m_name} in {time.time()-t_m:.1f}s", flush=True)
    return rows_m, cal_recs_m

def run_new_geology_phase():
    t_start = time.time()
    p1_spec = json.loads((SUITE_ROOT / 'core_p1/inputs/p1_experiment_spec.json').read_text())
    tasks = []
    for group_key, suite_prefix in [('development', 'DevMain'), ('heldout_new_geology', 'NewGeologyMain')]:
        for m_spec in p1_spec['models'][group_key]:
            tasks.append((group_key, suite_prefix, m_spec))

    rows_new, calibration_records = [], []
    with ProcessPoolExecutor(max_workers=2) as ex:
        results = list(ex.map(_run_single_new_model, tasks))
    for rm, cr in results:
        rows_new.extend(rm)
        calibration_records.extend(cr)

    res_dir = SUITE_ROOT / 'results'
    with (res_dir / 'metrics_p1_historical_paired.csv').open(encoding='utf-8-sig', newline='') as f:
        rows_hist = list(csv.DictReader(f))
    with (res_dir / 'metrics_p1_multical_extra.csv').open(encoding='utf-8-sig', newline='') as f:
        rows_mc = list(csv.DictReader(f))
    rows_all = rows_hist + rows_mc + rows_new
    fieldnames = list(rows_hist[0].keys())
    with (res_dir / 'metrics_p1_all_models_multical.csv').open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader(); w.writerows(rows_all)
    cal_hist = json.loads((res_dir / 'calibration_historical.json').read_text())
    (res_dir / 'calibration_estimates.json').write_text(json.dumps(cal_hist + calibration_records, ensure_ascii=False, indent=2), encoding='utf-8')

    eval_script = SUITE_ROOT / 'scripts/evaluate_candidate.py'
    for cand in ['CombinedUncertainty', 'NoiseUncertaintyOnly', 'QModelErrorOnly', 'StrongerDamping', 'CombinedUncertainty_NoStructure']:
        for suffix, csv_name in [('', 'metrics_p1_historical_paired.csv'), ('_all_models', 'metrics_p1_all_models_multical.csv')]:
            out_sub = res_dir / f'eval_{cand}_vs_Frozen{suffix}'
            if out_sub.exists(): shutil.rmtree(out_sub)
            subprocess.run([
                sys.executable, str(eval_script),
                '--metrics', str(res_dir / csv_name),
                '--candidate', cand, '--baseline', 'Frozen', '--route', 'Split',
                '--output', str(out_sub)
            ], check=True)
    print(f"New geology + evaluation phase completed in {time.time()-t_start:.1f}s!", flush=True)

if __name__ == '__main__':
    phase = sys.argv[1] if len(sys.argv) > 1 else 'historical'
    if phase == 'historical': run_historical_phase()
    elif phase == 'new_geology': run_new_geology_phase()
