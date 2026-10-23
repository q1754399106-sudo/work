"""Strict, standard-library experiment audit. CSV evidence never proves blindness.

Usage: python strict_evaluate.py --metrics metrics.csv --spec frozen_spec.json
       --metadata RUN_METADATA.json --output NEW_DIRECTORY

Absent metrics are reported as NOT_RUN/BLOCKED_ENGINE. Failed solve rows remain
in the paired matrix and cause FAIL; their NaN scores are never averaged away.
"""
from pathlib import Path
import argparse
import collections
import csv
import hashlib
import html
import json
import math
import random
import statistics

METRICS = ['RelativeError', 'EventAmplitudeError', 'WeakEventError',
           'FaultNeighborhoodError', 'PinchoutTailError', 'AVOInterceptError',
           'AVOGradientError']
KEYS = ['Model', 'Suite', 'SNRdB', 'CalibrationSeed', 'Seed', 'Scenario']
SUCCESS = {'COMPLETED', 'COMPLETE', 'SUCCESS', 'SUCCEEDED', 'OK', 'PASS'}
DEFAULT_POLICY = dict(main_equal_model_relative_error_improvement_min=.05,
    any_suite_model_mean_degradation_max=.02,
    event_and_weak_and_AVO_error_degradation_max=.05,
    local_near_zero_absolute_allowance=.001,
    pinchout_degradation_relative_max=.10,
    failure_fraction_error_gt_1_increase_max=0,
    normal_residual_max=2e-7)


class InvalidInput(ValueError):
    pass


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    try:
        return math.isfinite(float(value))
    except (ValueError, TypeError):
        return False


def require(condition, message):
    if not condition:
        raise InvalidInput(message)


def percentile(values, proportion):
    ordered = sorted(values)
    x = (len(ordered) - 1) * proportion
    lo = int(x)
    hi = min(lo + 1, len(ordered) - 1)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (x - lo)


def distribution(values):
    return dict(n=len(values), mean=statistics.mean(values),
        median=statistics.median(values), p90=percentile(values, .90),
        p95=percentile(values, .95), worst=max(values),
        failures=sum(x > 1 for x in values),
        failure_fraction=sum(x > 1 for x in values) / len(values))


def applicable(scenario):
    result = ['RelativeError', 'EventAmplitudeError']
    if scenario.endswith('AVO'):
        result += ['AVOInterceptError', 'AVOGradientError']
    if scenario == 'WeakSection':
        result.append('WeakEventError')
    if scenario == 'PinchoutSection':
        result.append('PinchoutTailError')
    if scenario == 'FaultSection':
        result.append('FaultNeighborhoodError')
    return result


def canonical_number(value):
    require(finite(value), 'Nonfinite SNR/key value: ' + str(value))
    return format(float(value), '.12g')


def is_failed(row):
    return row.get('Status', 'COMPLETED').upper() not in SUCCESS


def check_metadata(metadata, metadata_path, spec_path, spec, fixture):
    require(metadata.get('frozen_experiment_spec_sha256') == sha256(spec_path),
            'Frozen specification hash is absent or does not match --spec bytes')
    require(metadata.get('all_failures_retained') is True,
            'all_failures_retained must be true; expected matrix is checked separately')
    require(metadata.get('acceptance_thresholds_frozen_before_test') is True,
            'acceptance_thresholds_frozen_before_test must be true')
    require(bool(metadata.get('truth_isolation_audit')), 'Missing truth isolation audit record')
    require(bool(metadata.get('new_geological_holdout')), 'Missing geological holdout provenance')
    require(fixture or metadata.get('actual_matlab_execution') is True,
            'Nonempty scientific metrics require actual_matlab_execution=true')
    require(fixture or metadata.get('status') == 'COMPLETED',
            'Scientific run metadata must report COMPLETED; partial runs cannot pass')
    records = metadata.get('source_manifest')
    if records is None:
        hashes = metadata.get('source_hashes', {})
        paths = metadata.get('source_hash_paths', {})
        records = [dict(path=paths.get(k, k), sha256=v) for k, v in hashes.items()]
    require(isinstance(records, list) and records, 'Missing source manifest')
    seen = set()
    actual_records = {}
    for record in records:
        path = Path(record['path'])
        if not path.is_absolute():
            path = metadata_path.parent / path
        path = path.resolve()
        require(path not in seen, 'Duplicate source manifest path: ' + str(path))
        seen.add(path)
        require(path.is_file(), 'Missing recorded source: ' + str(path))
        require(sha256(path) == record['sha256'], 'Source hash mismatch: ' + str(path))
        actual_records[str(path)] = record['sha256']
    if not fixture:
        freeze = spec.get('freeze')
        require(isinstance(freeze, dict), 'Scientific test requires frozen source/config envelope')
        require(freeze.get('development_status') == 'COMPLETED'
                and freeze.get('development_actual_matlab_execution') is True,
                'Frozen specification requires completed actual development execution')
        snapshot = Path(freeze['config_snapshot'])
        if not snapshot.is_absolute():
            snapshot = spec_path.parent / snapshot
        require(snapshot.is_file() and sha256(snapshot) == freeze['config_sha256'],
                'Frozen configuration snapshot is absent or changed')
        original = json.loads(snapshot.read_text(encoding='utf-8-sig'))
        require(original == {k: v for k, v in spec.items() if k != 'freeze'},
                'Frozen spec differs from its original configuration snapshot')
        identity = freeze['config_sha256'] + '\n'
        frozen_records = {}
        for record in freeze['sources']:
            path = Path(record['path'])
            require(path.is_absolute() and path.is_file() and sha256(path) == record['sha256'],
                    'Frozen source changed or missing: ' + str(path))
            require(str(path.resolve()) not in frozen_records, 'Duplicate frozen source path')
            frozen_records[str(path.resolve())] = record['sha256']
            identity += record['path'] + '=' + record['sha256'] + '\n'
        require(hashlib.sha256(identity.encode('utf-8')).hexdigest() == freeze['expected_hash'],
                'Frozen source/config identity hash mismatch')
        require(frozen_records and all(frozen_records.get(path) == value
            for path, value in actual_records.items()),
            'Actual execution sources are not covered by frozen source manifest')
    return len(records)


def inspect_matrix(rows, spec):
    required = set(KEYS + METRICS + ['Method', 'Route', 'PriorHash',
        'SharedInformationHash', 'NormalResidual', 'FinalFlag'])
    require(rows and required.issubset(rows[0]),
            'Empty CSV or missing columns: ' + str(sorted(required - set(rows[0] if rows else {}))))
    test_models = list(spec.get('test_model_ids', []))
    development = list(spec.get('development_model_ids', []))
    require(len(set(test_models)) == len(test_models) and len(test_models) >= 4,
            'At least four unique test model IDs are required')
    require(len(set(development)) == len(development) and len(development) >= 2,
            'At least two unique development model IDs are required')
    require(not set(test_models) & set(development), 'Development/test model IDs overlap')
    scenarios = set(spec.get('scenario_names', []))
    require({'WeakSection', 'FaultSection', 'PinchoutSection'}.issubset(scenarios)
            and any(x.endswith('AVO') for x in scenarios), 'Required local scenarios absent from frozen spec')
    branches = list(spec.get('candidate_branches', []))
    require(len(branches) == len(set(branches)) and branches,
            'Frozen spec must list unique candidate_branches')
    controls = spec.get('controls', [])
    controls = [x['name'] if isinstance(x, dict) else x for x in controls]
    methods = ['Baseline', 'Frozen'] + branches + controls
    require(len(methods) == len(set(methods)), 'Overlapping method/control names')
    suites = spec.get('suite_snr_db')
    suite_configs = {item['name']: item for item in spec.get('suites', [])}
    if not suites and suite_configs:
        suites = {name: item['snr_db'] for name, item in suite_configs.items()}
    require(isinstance(suites, dict) and suites, 'Frozen spec needs suite_snr_db matrix')
    normalized_suites = {suite: {canonical_number(x) for x in (snrs if isinstance(snrs, list) else [snrs])}
                         for suite, snrs in suites.items()}
    require('Main' in suites, 'Main suite is required')
    require(set(spec.get('required_suites', suites)) == set(suites),
            'required_suites and suite_snr_db differ')
    require({canonical_number(x) for x in spec.get('snr_db', [])}.issubset(
        set.union(*normalized_suites.values())), 'Specified SNR missing from suite matrix')
    for suite in normalized_suites:
        config = suite_configs.setdefault(suite, dict(name=suite, scenario_names=sorted(scenarios),
            calibration_count=3, observation_count=6))
        if 'snr_db' in config:
            expected_snr = config['snr_db'] if isinstance(config['snr_db'], list) else [config['snr_db']]
            require({canonical_number(x) for x in expected_snr} == normalized_suites[suite],
                    'Suite SNR declarations disagree: ' + suite)
        require(set(config.get('scenario_names', [])).issubset(scenarios)
                and config.get('scenario_names'), 'Invalid frozen suite scenarios: ' + suite)
        if suite in {'Main', 'LowerSNR'}:
            require(config.get('calibration_count', 0) >= 3 and config.get('observation_count', 0) >= 6,
                    'Main/LowerSNR requires at least three calibrations and six observations')
            require({'SmoothAVO', 'WeakSection', 'FaultSection', 'PinchoutSection'}.issubset(
                    set(config['scenario_names'])),
                    'Main/LowerSNR must contain the four required white-noise scenarios')
    require(set.union(*(set(config['scenario_names']) for config in suite_configs.values())) == scenarios,
            'Frozen suite matrix does not cover all globally declared scenarios')
    lookup = {}
    failed_rows = []
    for source in rows:
        row = dict(source)
        row['SNRdB'] = canonical_number(row['SNRdB'])
        require(row['Model'] in test_models, 'Unexpected/non-holdout model: ' + row['Model'])
        require(row['Suite'] in normalized_suites and row['SNRdB'] in normalized_suites[row['Suite']],
                'Unexpected suite/SNR: ' + row['Suite'] + '/' + row['SNRdB'])
        require(row['Scenario'] in scenarios, 'Unexpected scenario: ' + row['Scenario'])
        require(row['Scenario'] in suite_configs[row['Suite']]['scenario_names'],
                'Scenario outside frozen suite plan: ' + row['Suite'] + '/' + row['Scenario'])
        require(row['Method'] in methods and row['Route'] in {'Direct', 'Split'},
                'Unexpected method/route: ' + row['Method'] + '/' + row['Route'])
        require(all(str(row[k]).strip() for k in KEYS + ['PriorHash', 'SharedInformationHash']),
                'Empty pairing key or provenance hash')
        key = tuple(row[k] for k in KEYS) + (row['Method'], row['Route'])
        require(key not in lookup, 'Duplicate case: ' + repr(key))
        failed = is_failed(row)
        if failed:
            failed_rows.append(dict(key=list(key), status=row.get('Status'),
                                    message=row.get('FailureMessage', '')))
        else:
            for metric in applicable(row['Scenario']) + ['NormalResidual', 'FinalFlag']:
                require(finite(row[metric]), 'Nonfinite successful metric ' + metric + ': ' + repr(key))
                require(float(row[metric]) >= 0, 'Negative error/solver metric ' + metric + ': ' + repr(key))
        lookup[key] = row
    observations = {key[:-2] for key in lookup}
    for observation in observations:
        information = set()
        for method in methods:
            pair = [observation + (method, route) for route in ['Direct', 'Split']]
            require(all(key in lookup for key in pair), 'Missing method/route at: ' + repr(observation + (method,)))
            direct, split = (lookup[key] for key in pair)
            require(direct['PriorHash'] == split['PriorHash'],
                    'Direct/Split prior hash differs: ' + repr(observation + (method,)))
            information.update(lookup[key]['SharedInformationHash'] for key in pair if not is_failed(lookup[key]))
        require(len(information) <= 1, 'Information differs across successful methods/routes: ' + repr(observation))
    # Main/LowerSNR repeat independent calibration; frozen one-factor stress
    # suites may deliberately use fewer calibration layers and local scenes.
    calibration_seeds = {}
    observation_seeds = {}
    planned_models = {item['id']: item for item in spec.get('models', [])}
    for model in test_models:
        modelkeys = [key for key in observations if key[0] == model]
        calibrations = {key[3] for key in modelkeys if key[1] in {'Main', 'LowerSNR'}}
        require(len(calibrations) >= max(3, spec.get('calibration_realizations_per_model_min', 3)),
                'Insufficient independent calibrations: ' + model)
        calibration_seeds[model] = sorted(calibrations)
        plan = planned_models.get(model)
        if planned_models:
            require(plan is not None, 'Test model absent from frozen model table: ' + model)
            require(calibrations == {str(x) for x in plan['calibration_seeds']},
                    'Calibration seeds differ from frozen model plan: ' + model)
        for calibration in calibrations:
            seeds = {key[4] for key in modelkeys if key[3] == calibration and key[1] in {'Main', 'LowerSNR'}}
            require(len(seeds) >= max(6, spec.get('observation_seeds_per_calibration_min', 6)),
                    'Insufficient observation seeds: ' + model + '/' + calibration)
            observation_seeds[(model, calibration)] = sorted(seeds)
            require(calibration not in seeds, 'Calibration and observation seed are identical')
        for suite, snrs in normalized_suites.items():
            config = suite_configs[suite]
            suitekeys = [key for key in modelkeys if key[1] == suite]
            suite_calibrations = {key[3] for key in suitekeys}
            require(len(suite_calibrations) == config['calibration_count'],
                    'Missing/extra calibration layer in frozen suite: ' + model + '/' + suite)
            if suite in {'Main', 'LowerSNR'}:
                require(suite_calibrations == calibrations, 'Main/LowerSNR calibration layers differ: ' + model)
            else:
                require(suite_calibrations.issubset(calibrations), 'Stress uses unknown calibration realization')
            for calibration in suite_calibrations:
                seeds = {key[4] for key in suitekeys if key[3] == calibration}
                require(len(seeds) == config['observation_count'],
                        'Missing/extra observation layer in frozen suite: ' + model + '/' + suite + '/' + calibration)
                if suite in {'Main', 'LowerSNR'}:
                    require(set(observation_seeds[(model, calibration)]) == seeds,
                            'Main/LowerSNR observation seed layers differ')
                if plan:
                    planned_cals = [str(x) for x in plan['calibration_seeds']]
                    require(suite_calibrations == set(planned_cals[:config['calibration_count']]),
                            'Suite calibration seeds differ from frozen order: ' + model + '/' + suite)
                    index = planned_cals.index(calibration)
                    expected_seeds = {str(x) for x in plan['observation_seeds'][index][:config['observation_count']]}
                    require(seeds == expected_seeds, 'Observation seeds differ from frozen model plan: ' + model)
                for snr in snrs:
                    for seed in seeds:
                        for scenario in config['scenario_names']:
                            require((model, suite, snr, calibration, seed, scenario) in observations,
                                'Missing model/SNR/calibration/seed/scenario case: ' +
                                repr((model, suite, snr, calibration, seed, scenario)))
    return lookup, observations, methods, branches, failed_rows, calibration_seeds, observation_seeds


def audit_declared_seeds(metadata, calibration_seeds, observation_seeds):
    declared = metadata.get('calibration_seeds_by_model')
    if declared:
        for model, seeds in calibration_seeds.items():
            require(model in declared and set(map(str, declared[model])) == set(seeds),
                    'Declared calibration seeds differ: ' + model)
    declared = metadata.get('observation_seeds_by_calibration')
    if declared:
        for (model, calibration), seeds in observation_seeds.items():
            # Accept model->{calibration: seeds} or model/calibration->seeds.
            entry = declared.get(model + '/' + calibration)
            if entry is None and isinstance(declared.get(model), dict):
                entry = declared[model].get(calibration, declared[model].get('cal_' + calibration))
            require(entry is not None and set(map(str, entry)) == set(seeds),
                    'Declared observation seeds differ: ' + model + '/' + calibration)


def bootstrap(pairs, spec):
    # A seed block averages scenes; geologies, then calibrations, then observation
    # seed blocks are resampled. Pairing is preserved at every level.
    tree = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    for baseline, candidate in pairs:
        tree[baseline['Model']][baseline['CalibrationSeed']][baseline['Seed']].append(
            (float(baseline['RelativeError']), float(candidate['RelativeError'])))
    rng = random.Random(int(spec.get('bootstrap_seed', 20261019)))
    count = int(spec.get('bootstrap_replicates', 1000))
    require(100 <= count <= 20000, 'bootstrap_replicates must lie in [100,20000]')
    geologies = sorted(tree)
    gains = []
    differences = []
    for _ in range(count):
        geological_means = []
        for model in rng.choices(geologies, k=len(geologies)):
            calibrations = sorted(tree[model])
            calibration_means = []
            for calibration in rng.choices(calibrations, k=len(calibrations)):
                seeds = sorted(tree[model][calibration])
                blocks = [tree[model][calibration][seed] for seed in rng.choices(seeds, k=len(seeds))]
                calibration_means.append(tuple(statistics.mean(
                    statistics.mean(x[i] for x in block) for block in blocks) for i in [0, 1]))
            geological_means.append(tuple(statistics.mean(x[i] for x in calibration_means) for i in [0, 1]))
        b, c = (statistics.mean(x[i] for x in geological_means) for i in [0, 1])
        differences.append(c - b)
        if b > 0:
            gains.append(1 - c / b)
    return dict(kind='descriptive hierarchical paired bootstrap; geology/calibration/observation blocks',
        confidence_level=.95, replicates=count, geological_models=len(geologies),
        relative_improvement_interval=[percentile(gains, .025), percentile(gains, .975)] if gains else None,
        error_difference_interval=[percentile(differences, .025), percentile(differences, .975)],
        warning='Small geological sample; not a significance test or proof of generalization.')


def compare(lookup, observations, method, route, policy, spec):
    pairs = [(lookup[key + ('Frozen', route)], lookup[key + (method, route)]) for key in sorted(observations)]
    checks = []
    def check(name, passed, detail):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))
    failed = [(b, c) for b, c in pairs if is_failed(b) or is_failed(c)]
    solver = [row for pair in pairs for row in pair]
    check('solver_quality_and_all_failures_retained', all(not is_failed(row)
        and finite(row['NormalResidual']) and finite(row['FinalFlag'])
        and float(row['FinalFlag']) == 0
        and 0 <= float(row['NormalResidual']) <= min(2e-7, policy['normal_residual_max']) for row in solver),
        dict(failed_pair_count=len(failed), rows=len(solver)))
    if failed:
        return dict(status='FAIL', method=method, route=route, checks=checks,
                    unavailable_metrics='Failed solve rows retained; summaries unavailable to prevent selective deletion.',
                    failed_pairs=len(failed))
    groups = []
    groupkeys = sorted({(b['Suite'], b['Model'], b['SNRdB']) for b, _ in pairs})
    for suite, model, snr in groupkeys:
        subset = [(b, c) for b, c in pairs if (b['Suite'], b['Model'], b['SNRdB']) == (suite, model, snr)]
        bstats, cstats = (distribution([float(pair[i]['RelativeError']) for pair in subset]) for i in [0, 1])
        groups.append(dict(Suite=suite, Model=model, SNRdB=snr, baseline=bstats, candidate=cstats))
        prefix = suite + '/' + model + '/SNR' + snr
        check(prefix + '/mean_noninferiority', cstats['mean'] <= bstats['mean'] *
              (1 + policy['any_suite_model_mean_degradation_max']), groups[-1])
        check(prefix + '/failure_fraction', cstats['failure_fraction'] <= bstats['failure_fraction'] +
              policy['failure_fraction_error_gt_1_increase_max'],
              dict(baseline=bstats['failure_fraction'], candidate=cstats['failure_fraction']))
        for scenario in sorted({b['Scenario'] for b, _ in subset}):
            local = [(b, c) for b, c in subset if b['Scenario'] == scenario]
            for metric in applicable(scenario)[1:]:
                x, y = (statistics.mean(float(pair[i][metric]) for pair in local) for i in [0, 1])
                rel = policy['pinchout_degradation_relative_max'] if metric == 'PinchoutTailError' else policy['event_and_weak_and_AVO_error_degradation_max']
                allowance = max(x * rel, policy['local_near_zero_absolute_allowance'])
                check(prefix + '/' + scenario + '/' + metric, y <= x + allowance,
                      dict(baseline=x, candidate=y, allowance=allowance, n=len(local)))
    main = [(b, c) for b, c in pairs if b['Suite'] == 'Main']
    require(main, 'No Main data')
    mainmeans = [tuple(statistics.mean(float(pair[i]['RelativeError']) for pair in main
        if pair[0]['Model'] == model) for i in [0, 1]) for model in spec['test_model_ids']]
    b, c = (statistics.mean(x[i] for x in mainmeans) for i in [0, 1])
    gain = 1 - c / b if b > 0 else None
    check('main_equal_model_improvement', gain is not None and gain >= policy['main_equal_model_relative_error_improvement_min'],
          dict(baseline=b, candidate=c, relative_improvement=gain))
    diagnostics = collections.Counter()
    for _, row in pairs:
        if finite(row.get('ScaleEndpoint')) and float(row['ScaleEndpoint']) != 0:
            diagnostics['scale_endpoints'] += 1
        if finite(row.get('TargetMismatch')) and float(row['TargetMismatch']) > .005:
            diagnostics['target_mismatches'] += 1
    return dict(status='PROVISIONAL_PASS' if all(x['passed'] for x in checks) else 'FAIL',
        method=method, route=route, pairs=len(pairs), equal_model_main_improvement=gain,
        groups=groups, checks=checks, diagnostics=dict(diagnostics),
        hierarchical_uncertainty=bootstrap(main, spec))


def write_reports(output, result):
    (output / 'acceptance.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    lines = ['# P1 补偿实验：严格评价', '', '**结论：' + result['status'] + '。**', '',
        'CSV与元数据校验不证明盲测、真值隔离或地质独立性；数值通过仅为暂定工程结论，不自动替换Frozen。']
    if result.get('artificial_fixture'):
        lines += ['', '**人工测试夹具：不是反演结果或科学证据。**']
    if result.get('reason'):
        lines += ['', result['reason']]
    for method, branch in result.get('branches', {}).items():
        lines += ['', '## ' + method, '', '主路线 Split：' + branch['status'], '',
            '| Route / 路线 | Status / 结论 | Main gain / 主试验改善 |', '|---|---|---:|']
        for route, comparison in branch['routes'].items():
            gain = comparison.get('equal_model_main_improvement')
            lines += ['|' + route + '|' + comparison['status'] + '|' + (f'{gain:.2%}' if gain is not None else '不可计算') + '|']
        for check in branch['routes']['Split']['checks']:
            if not check['passed']:
                lines += ['', '- 未通过：' + check['name'] + '；' + json.dumps(check['detail'], ensure_ascii=False)]
    lines += ['', '分层区间为描述性结果；按地质、校准与观测块重采样。完整分布、局部指标、数值检查及失败记录见acceptance.json。']
    markdown = '\n'.join(lines) + '\n'
    (output / 'EFFECTS_ZH.md').write_text(markdown, encoding='utf-8')
    body = '<h1>P1 effects / P1效果</h1><pre>' + html.escape(markdown) + '</pre>'
    body += '<h2>Complete audit / 完整审计</h2><pre>' + html.escape(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)) + '</pre>'
    (output / 'effects.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>P1 effects</title><style>body{max-width:1200px;margin:32px auto;padding:20px;font:16px/1.6 system-ui}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style>' + body + '</html>', encoding='utf-8')


def evaluate(metrics, spec_path, metadata_path, output):
    metrics, spec_path, metadata_path, output = map(Path, (metrics, spec_path, metadata_path, output))
    require(not output.exists(), 'Use a fresh output directory')
    output.mkdir(parents=True)
    result = dict(status='INVALID_INPUT', scope='Strict numerical/provenance audit; external scientific audit remains required.')
    try:
        spec = json.loads(spec_path.read_text(encoding='utf-8-sig'))
        metadata = json.loads(metadata_path.read_text(encoding='utf-8-sig'))
        fixture = metadata.get('artificial_fixture') is True
        result['artificial_fixture'] = fixture
        if not metrics.is_file() or metrics.stat().st_size == 0:
            status = 'BLOCKED_ENGINE' if 'BLOCKED_ENGINE' in str(metadata.get('status', '')).upper() else 'NOT_RUN'
            result.update(status=status, reason='没有实际实验CSV；未生成或推断实验效果。',
                          actual_matlab_execution=metadata.get('actual_matlab_execution', False))
            write_reports(output, result)
            return result
        with metrics.open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            result.update(status='NOT_RUN', reason='CSV只有表头，没有实际观测；未生成实验效果。')
            write_reports(output, result)
            return result
        sources = check_metadata(metadata, metadata_path, spec_path, spec, fixture)
        lookup, observations, methods, branches, failures, calseeds, obsseeds = inspect_matrix(rows, spec)
        audit_declared_seeds(metadata, calseeds, obsseeds)
        policy = dict(DEFAULT_POLICY)
        policy.update(spec.get('acceptance', {}))
        comparisons = {}
        for method in branches:
            routes = {route: compare(lookup, observations, method, route, policy, spec) for route in ['Direct', 'Split']}
            direct_numerical_ok = routes['Direct']['checks'][0]['passed']
            comparisons[method] = dict(status=routes['Split']['status'] if direct_numerical_ok else 'FAIL',
                                      matched_prior_direct_numerical_check=direct_numerical_ok, routes=routes)
        # Controls are scored identically and reported without changing branch eligibility.
        controls = {method: {route: compare(lookup, observations, method, route, policy, spec)
            for route in ['Direct', 'Split']} for method in methods if method not in branches + ['Baseline', 'Frozen']}
        baseline_ok = all(not is_failed(row) and finite(row['NormalResidual'])
            and finite(row['FinalFlag']) and float(row['FinalFlag']) == 0
            and 0 <= float(row['NormalResidual']) <= min(2e-7, policy['normal_residual_max'])
            for key, row in lookup.items() if key[-2] == 'Baseline')
        passing = [name for name, branch in comparisons.items() if branch['status'] == 'PROVISIONAL_PASS']
        failing = [name for name in comparisons if name not in passing]
        if not baseline_ok or not passing:
            status = 'FAIL'
        elif failing:
            status = 'ARTIFICIAL_MIXED_RESULTS' if fixture else 'MIXED_RESULTS'
        else:
            status = 'ARTIFICIAL_FIXTURE_PASS' if fixture else 'PROVISIONAL_PASS'
        result.update(status=status, passing_branches=passing, failing_branches=failing,
            rows=len(rows), observations=len(observations), verified_source_files=sources,
            metrics_sha256=sha256(metrics), frozen_experiment_spec_sha256=sha256(spec_path),
            metadata_sha256=sha256(metadata_path), branches=comparisons, controls=controls,
            baseline_control_numerically_valid=baseline_ok,
            failed_rows=failures, source_and_spec_hashes_verified=True,
            external_requirements=['Independent geological holdout audit', 'Truth-isolation audit of inverse-only inputs',
                'Pre-test freeze chronology and tuning-budget audit', 'Physical calibration independence, not merely distinct seed labels'],
            metadata_claims=dict(truth_isolation_audit=metadata.get('truth_isolation_audit'),
                new_geological_holdout=metadata.get('new_geological_holdout')),
            decision='Do not automatically replace Frozen. Review each branch, Direct comparator, controls and local costs.')
        paired = []
        for method in branches:
            for route in ['Direct', 'Split']:
                for key in sorted(observations):
                    baseline, candidate = lookup[key + ('Frozen', route)], lookup[key + (method, route)]
                    record = dict(zip(KEYS, key), Method=method, Route=route,
                                  BaselineStatus=baseline.get('Status', 'COMPLETED'), CandidateStatus=candidate.get('Status', 'COMPLETED'))
                    for metric in METRICS:
                        b, c = baseline[metric], candidate[metric]
                        record.update({metric + '_baseline': b, metric + '_candidate': c,
                                       metric + '_difference': float(c) - float(b) if finite(b) and finite(c) else ''})
                    paired.append(record)
        with (output / 'paired_differences.csv').open('w', encoding='utf-8-sig', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(paired[0]))
            writer.writeheader()
            writer.writerows(paired)
    except (InvalidInput, OSError, ValueError, KeyError, TypeError) as error:
        result.update(status='INVALID_INPUT', reason=str(error))
    write_reports(output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['metrics', 'spec', 'metadata', 'output']:
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    try:
        result = evaluate(args.metrics, args.spec, args.metadata, args.output)
    except InvalidInput as error:
        parser.error(str(error))
    print(json.dumps({k: result[k] for k in ['status', 'reason'] if k in result}, ensure_ascii=False))
    return 2 if result['status'] == 'INVALID_INPUT' else 0


if __name__ == '__main__':
    raise SystemExit(main())
