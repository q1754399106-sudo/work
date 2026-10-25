"""Post-hoc diagnostics for a COMPLETED test; never changes raw run evidence.

python audit_solver_logs.py --run experiments/test_v2/matlab --output NEW_DIR
python audit_solver_logs.py --self-test

This is a reporting supplement, not an acceptance evaluator or tuning tool.
In particular the original CSV TargetMismatch is the selected Direct mismatch
on both routes. The selected Split mismatch is reconstructed from its own log.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import tempfile
import unittest

CASE_KEYS = ('Model', 'Suite', 'Scenario', 'CalibrationSeed', 'Seed', 'Method')
NUMERIC_RECORD_FIELDS = ('flag', 'initial_flag', 'iterations', 'normal_residual',
                         'data_residual')
RECORD_FIELDS = NUMERIC_RECORD_FIELDS + ('pcg_relative_residual', 'scale',
    'mahalanobis_quadratic', 'fast_normal_residual_before_restart')
LEGACY_METHODS = {'Baseline', 'Frozen'}


def number(value):
    if value is None or isinstance(value, (list, dict)):
        return None
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    return value


def as_list(value):
    return value if isinstance(value, list) else ([value] if isinstance(value, dict) else [])


def key_of(value):
    return tuple(str(value.get(k, '')) for k in CASE_KEYS)


def identity(key):
    return dict(zip(CASE_KEYS, key))


def mismatch(residual, target):
    r, t = number(residual), number(target)
    return abs(r / t - 1) if r is not None and t is not None and t > 0 else None


def selected_target(log, route):
    if log.get('Method') in LEGACY_METHODS:
        # solve_v2 reports route-specific legacy Euclidean residuals. Both are
        # compared to inherited log.target, not a new whitened target.
        return number(log.get('target')), 'legacy Euclidean residual / inherited log.target'
    d = log.get('diagnostics', {})
    squared = number(d.get(route.lower() + '_target_squared')) if isinstance(d, dict) else None
    return (math.sqrt(squared) if squared is not None and squared > 0 else None,
            'fixed-covariance Mahalanobis residual / route-specific sqrt(target_squared)')


def record_diagnostics(record, threshold):
    if not isinstance(record, dict):
        record = {}
    bad = [name for name in NUMERIC_RECORD_FIELDS if number(record.get(name)) is None]
    # Preserve any optional null/nonfinite numeric diagnostics as well.
    bad += [name for name in RECORD_FIELDS if name in record
            and name not in bad and number(record[name]) is None]
    initial, flag = number(record.get('initial_flag')), number(record.get('flag'))
    normal = number(record.get('normal_residual'))
    solver = str(record.get('solver', ''))
    fallback = 'restart' in solver.lower() or 'cholesky' in solver.lower()
    trigger = initial is not None and initial != 0
    before = number(record.get('fast_normal_residual_before_restart'))
    return dict(NonfiniteOrMissingFields=';'.join(bad),
        FinalFlagNonzero=flag is not None and flag != 0,
        InitialFlagNonzero=trigger,
        NormalResidualExceeded=normal is not None and (normal < 0 or normal > threshold),
        RestartOrFallbackRecorded=fallback,
        InitialOrFastResidualRestartTrigger=trigger or (before is not None and before > threshold),
        NegativePrecisionQuadratic=record.get('negative_precision_quadratic') is True)


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    if not fields:
        fields = ['NoRecords']
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: ('' if v is None else v) for k, v in clean(row).items()} for row in rows)


def summarize_records(rows):
    normals = [number(r.get('normal_residual')) for r in rows]
    normals = [v for v in normals if v is not None]
    return dict(records_total=len(rows),
        nonzero_final_flags=sum(bool(r['FinalFlagNonzero']) for r in rows),
        nonzero_initial_flags=sum(bool(r['InitialFlagNonzero']) for r in rows),
        normal_residual_exceedances=sum(bool(r['NormalResidualExceeded']) for r in rows),
        nonfinite_or_missing_records=sum(bool(r['NonfiniteOrMissingFields']) for r in rows),
        restart_or_fallback_records=sum(bool(r['RestartOrFallbackRecorded']) for r in rows),
        restart_trigger_records=sum(bool(r['InitialOrFastResidualRestartTrigger']) for r in rows),
        negative_precision_quadratics=sum(bool(r['NegativePrecisionQuadratic']) for r in rows),
        max_normal_residual=max(normals) if normals else None,
        solver_counts=dict(collections.Counter(str(r.get('solver', 'MISSING')) for r in rows)))


def run_audit(run, output):
    run, output = Path(run).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError('Output must be a fresh directory')
    if output == run or run in output.parents:
        raise ValueError('Write the post-hoc audit outside the raw run directory')
    manifest, issues = [], []

    def source(path):
        content = path.read_bytes()
        manifest.append(dict(Path=str(path), SHA256=hashlib.sha256(content).hexdigest(), Bytes=len(content)))
        return content

    def read_json(path):
        value = json.loads(source(path).decode('utf-8-sig'))
        if not isinstance(value, dict):
            raise ValueError('Expected JSON object: ' + str(path))
        return value

    # Reject partial results before reading any solver/calibration evidence or
    # creating an output directory. No interim evaluation is performed.
    metadata = read_json(run / 'RUN_METADATA.json')
    status = read_json(run / 'run_status.json')
    if metadata.get('status') != 'COMPLETED' or status.get('status') != 'COMPLETED':
        raise ValueError('Audit refused: test has not COMPLETED')
    if metadata.get('stage') != 'test' or status.get('stage') != 'test':
        raise ValueError('Audit refused: only completed held-out test runs are accepted')
    if metadata.get('actual_matlab_execution') is not True:
        raise ValueError('Actual MATLAB execution is required')
    spec = read_json(run / 'experiment_spec_used.json')
    threshold = min(2e-7, float(spec.get('acceptance', {}).get('normal_residual_max', 2e-7)))
    with (run / 'metrics.csv').open(encoding='utf-8-sig', newline='') as stream:
        metrics = list(csv.DictReader(stream))
    source(run / 'metrics.csv')
    if not metrics or len(metrics) != metadata.get('rows') or len(metrics) != status.get('rows'):
        raise ValueError('Completed metrics row count does not match run metadata/status')
    expected = collections.defaultdict(list)
    for row in metrics:
        expected[key_of(row)].append(row)
    expected_names = {f"{k[0]}_{k[1]}_{k[2]}_cal{k[3]}_seed{k[4]}_{k[5]}.json": k for k in expected}
    records, selected, cases, seen = [], [], [], collections.Counter()
    corrupt_solver = 0

    def case_rows(log, key, path, state='PARSED'):
        common = dict(identity(key), LogPath=str(path), LogStatus=state)
        case_issues = []
        matched = expected.get(key, [])
        if not matched:
            case_issues.append('UnexpectedSolverCase')
        by_route = {route: [r for r in matched if r.get('Route') == route] for route in ('Direct', 'Split')}
        if any(len(v) != 1 for v in by_route.values()):
            case_issues.append('MissingOrDuplicateMetricRoute')
        history = as_list(log.get('records'))
        grids = log.get('scale_grid', [])
        grids = grids if isinstance(grids, list) else [grids]
        if len(history) != len(grids) + 1:
            case_issues.append('HistoryCountNotSearchGridPlusOneSplit')
        if log.get('FailureMessage'):
            case_issues.append('LoggedFailureMessage')
        for index, record in enumerate(history, 1):
            record = record if isinstance(record, dict) else {}
            purpose = 'SearchedDirect' if index <= len(grids) else ('FinalSplit' if index == len(grids) + 1 else 'UnexpectedExtraRecord')
            rec = dict(common, RecordIndex=index, Purpose=purpose,
                SearchGridScale=grids[index-1] if index <= len(grids) else None,
                **{k: record.get(k) for k in RECORD_FIELDS}, solver=record.get('solver', ''),
                **record_diagnostics(record, threshold), RawRecordJSON=json.dumps(clean(record), ensure_ascii=False))
            records.append(rec)
        direct_mismatch, split_mismatch = None, None
        route_rows = []
        for route in ('Direct', 'Split'):
            record = log.get(route.lower() + '_record', {})
            record = record if isinstance(record, dict) else {}
            if not record:
                case_issues.append('MissingSelected' + route + 'Record')
            target, units = selected_target(log, route)
            gap = mismatch(record.get('data_residual'), target)
            if gap is None:
                case_issues.append('Unavailable' + route + 'TargetMismatch')
            if route == 'Direct':
                direct_mismatch = gap
            else:
                split_mismatch = gap
            metric = by_route[route][0] if len(by_route[route]) == 1 else {}
            route_row = dict(common, Route=route, MetricStatus=metric.get('Status', 'MISSING'),
                Scale=log.get('scale'), ScaleEndpoint=log.get('endpoint'),
                Target=target, TargetUnits=units, ActualRouteTargetMismatch=gap,
                TargetMismatchGreaterThan0_005=gap is not None and gap > .005,
                OriginalCSVTargetMismatch=metric.get('TargetMismatch'),
                OriginalCSVTargetMismatchMeaning='selected Direct mismatch copied to both routes',
                OriginalLogDirectMismatch=log.get('residual_mismatch'),
                **{k: record.get(k) for k in RECORD_FIELDS}, solver=record.get('solver', ''),
                **record_diagnostics(record, threshold))
            selected.append(route_row)
            route_rows.append(route_row)
        if direct_mismatch is not None and number(log.get('residual_mismatch')) is not None:
            if not math.isclose(direct_mismatch, float(log['residual_mismatch']), abs_tol=1e-10, rel_tol=1e-8):
                case_issues.append('RecomputedDirectMismatchDiffersFromLog')
        summary = summarize_records(records[-len(history):]) if history else summarize_records([])
        cases.append(dict(common, FailureMessage=log.get('FailureMessage', ''),
            MatchedMetricRows=len(matched), SearchGridCount=len(grids), HistoryRecordCount=len(history),
            SelectedScale=log.get('scale'), ScaleEndpoint=log.get('endpoint'),
            SelectedDirectTargetMismatch=direct_mismatch, SelectedSplitTargetMismatch=split_mismatch,
            InitialFlagsNonzero=summary['nonzero_initial_flags'], FinalFlagsNonzero=summary['nonzero_final_flags'],
            NormalResidualExceeded=summary['normal_residual_exceedances'],
            NonfiniteOrMissingRecords=summary['nonfinite_or_missing_records'],
            RestartOrFallbackRecords=summary['restart_or_fallback_records'],
            Issues=';'.join(case_issues)))
        if case_issues:
            issues.append(dict(kind='SOLVER_CASE', **identity(key), path=str(path), details=case_issues))

    for path in sorted((run / 'solver_logs').glob('*.json')):
        try:
            log = read_json(path)
            key = key_of(log)
            seen[key] += 1
            case_rows(log, key, path)
        except (OSError, ValueError, TypeError) as exc:
            corrupt_solver += 1
            key = expected_names.get(path.name, ('', '', '', '', '', ''))
            seen[key] += 1
            case_rows({'FailureMessage': str(exc)}, key, path, 'CORRUPT')
    missing = [key for key in expected if key not in seen]
    for key in missing:
        case_rows({'FailureMessage': 'Expected solver log missing'}, key,
                  run / 'solver_logs' / next(name for name, k in expected_names.items() if k == key), 'MISSING')
    duplicates = [dict(identity(k), count=v) for k, v in seen.items() if v > 1]
    if duplicates:
        issues.append(dict(kind='DUPLICATE_SOLVER_LOG_CASES', cases=duplicates))

    calibration_rows, calibration_files, corrupt_calibration = [], [], 0
    expected_calibrations = {(r['Model'], r['Suite'], r['CalibrationSeed']) for r in metrics}
    expected_bootstrap = int(spec.get('algorithm', {}).get('bootstrap_replicates', 8))
    seen_calibrations = collections.Counter()

    def calibration_case(data, path, file_state='PARSED'):
        key = tuple(str(data.get(k, '')) for k in ('Model', 'Suite', 'CalibrationSeed'))
        common = dict(zip(('Model', 'Suite', 'CalibrationSeed'), key), LogPath=str(path), LogStatus=file_state,
                      BootstrapSeed=data.get('bootstrap_seed'))
        seen_calibrations[key] += 1
        if key not in expected_calibrations:
            issues.append(dict(kind='UNEXPECTED_CALIBRATION_CASE', **common))
        draws = as_list(data.get('bootstrap_estimates'))
        if len(draws) != expected_bootstrap:
            issues.append(dict(kind='CALIBRATION_BOOTSTRAP_COUNT', **common,
                               expected=expected_bootstrap, actual=len(draws)))
        estimates = [('Reference', 0, data.get('estimate'))] + [('Bootstrap', i+1,
            draws[i] if i < len(draws) else None) for i in range(max(expected_bootstrap, len(draws)))]
        for role, index, est in estimates:
            est = est if isinstance(est, dict) else {}
            fields = ('effective_Q', 'estimated_frequency', 'nominal_frequency', 'minimum_loss',
                'optimizer_exitflag', 'optimizer_evaluations', 'q_boundary', 'frequency_boundary')
            bad = [k for k in ('effective_Q', 'estimated_frequency', 'minimum_loss',
                              'optimizer_exitflag', 'optimizer_evaluations') if number(est.get(k)) is None]
            if not all(isinstance(est.get(k), bool) for k in ('q_boundary', 'frequency_boundary')):
                bad.append('boundary_flags')
            exitflag = number(est.get('optimizer_exitflag'))
            calibration_rows.append(dict(common, Role=role, BootstrapIndex=index,
                EstimateStatus='PARSED' if est else 'MISSING',
                **{k: est.get(k) for k in fields}, NonfiniteOrMissingFields=';'.join(bad),
                OptimizerNonpositiveExitFlag=exitflag is not None and exitflag <= 0,
                FailureMessage=data.get('FailureMessage', ''),
                EstimateSHA256=hashlib.sha256(json.dumps(clean(est), sort_keys=True).encode()).hexdigest() if est else ''))
        calibration_files.append(dict(common, BootstrapEstimatesPresent=len(draws),
            TruthWasUsedForEstimation=data.get('truth_was_used_for_estimation')))

    for path in sorted((run / 'calibration').glob('*.json')):
        try:
            calibration_case(read_json(path), path)
        except (OSError, ValueError, TypeError) as exc:
            corrupt_calibration += 1
            match = next((k for k in expected_calibrations if path.name == f'{k[0]}_{k[1]}_cal{k[2]}.json'), ('', '', ''))
            data = dict(zip(('Model', 'Suite', 'CalibrationSeed'), match), FailureMessage=str(exc))
            calibration_case(data, path, 'CORRUPT')
    missing_calibrations = expected_calibrations - set(seen_calibrations)
    for key in sorted(missing_calibrations):
        data = dict(zip(('Model', 'Suite', 'CalibrationSeed'), key), FailureMessage='Expected calibration log missing')
        calibration_case(data, run / 'calibration' / f'{key[0]}_{key[1]}_cal{key[2]}.json', 'MISSING')
    selected_summary = []
    for method, route in sorted({(r['Method'], r['Route']) for r in selected}):
        rows = [r for r in selected if r['Method'] == method and r['Route'] == route]
        gaps = [r['ActualRouteTargetMismatch'] for r in rows if r['ActualRouteTargetMismatch'] is not None]
        selected_summary.append(dict(Method=method, Route=route, cases=len(rows),
            nonzero_final_flags=sum(r['FinalFlagNonzero'] for r in rows),
            normal_residual_exceedances=sum(r['NormalResidualExceeded'] for r in rows),
            nonfinite_or_missing_records=sum(bool(r['NonfiniteOrMissingFields']) for r in rows),
            scale_endpoints=sum(r['ScaleEndpoint'] is True for r in rows),
            unavailable_target_mismatches=len(rows)-len(gaps),
            target_mismatches_gt_0_005=sum(v > .005 for v in gaps),
            mean_target_mismatch=statistics.mean(gaps) if gaps else None,
            max_target_mismatch=max(gaps) if gaps else None,
            target_units=sorted({r['TargetUnits'] for r in rows})))
    calibration_summary = dict(files=len(calibration_files),
        reference_estimates=sum(r['Role'] == 'Reference' for r in calibration_rows),
        bootstrap_estimates=sum(r['Role'] == 'Bootstrap' for r in calibration_rows),
        reference_estimates_present=sum(r['Role'] == 'Reference' and r['EstimateStatus'] == 'PARSED' for r in calibration_rows),
        bootstrap_estimates_present=sum(r['Role'] == 'Bootstrap' and r['EstimateStatus'] == 'PARSED' for r in calibration_rows),
        unique_model_calibration_seed_blocks=len({(r['Model'], r['CalibrationSeed']) for r in calibration_files}),
        unique_model_bootstrap_seed_blocks=len({(r['Model'], r['BootstrapSeed']) for r in calibration_files}),
        geological_models=len({r['Model'] for r in calibration_files if r['Model']}),
        interpretation='File and bootstrap counts are computational records; seed blocks repeat across suites. They are not independent geologies or independent calibration sample counts.')
    for role in ('Reference', 'Bootstrap'):
        rows = [r for r in calibration_rows if r['Role'] == role]
        prefix = role.lower() + '_'
        calibration_summary.update({prefix+'q_boundaries': sum(r['q_boundary'] is True for r in rows),
            prefix+'frequency_boundaries': sum(r['frequency_boundary'] is True for r in rows),
            prefix+'optimizer_nonpositive_exitflags': sum(r['OptimizerNonpositiveExitFlag'] for r in rows),
            prefix+'nonfinite_or_missing_estimates': sum(bool(r['NonfiniteOrMissingFields']) for r in rows),
            prefix+'optimizer_exitflag_counts': dict(collections.Counter(str(r['optimizer_exitflag']) for r in rows))})
    solver_summary = summarize_records(records)
    counts = dict(metrics_rows=len(metrics), expected_solver_logs=len(expected), solver_log_cases=len(cases),
        missing_solver_logs=len(missing), corrupt_solver_logs=corrupt_solver,
        duplicate_solver_log_cases=len(duplicates), expected_calibration_logs=len(expected_calibrations),
        calibration_logs=len(calibration_files), missing_calibration_logs=len(missing_calibrations),
        corrupt_calibration_logs=corrupt_calibration,
        duplicate_calibration_cases=sum(v > 1 for v in seen_calibrations.values()))
    findings = bool(issues) or any(solver_summary[k] for k in ('nonzero_final_flags',
        'nonzero_initial_flags', 'restart_or_fallback_records', 'nonfinite_or_missing_records',
        'normal_residual_exceedances', 'negative_precision_quadratics'))
    findings = findings or any(calibration_summary[k] for k in calibration_summary
                              if k.endswith(('boundaries', 'nonpositive_exitflags', 'missing_estimates')))
    findings = findings or any(r['target_mismatches_gt_0_005'] or r['scale_endpoints']
                              or r['unavailable_target_mismatches'] for r in selected_summary)
    findings = findings or any(counts[k] for k in ('missing_solver_logs', 'corrupt_solver_logs',
        'duplicate_solver_log_cases', 'missing_calibration_logs', 'corrupt_calibration_logs',
        'duplicate_calibration_cases'))
    input_hashes = {r['Path']: r['SHA256'] for r in manifest}
    result = dict(schema_version=1, status='COMPLETED_WITH_FINDINGS' if findings else 'COMPLETED_POSTHOC_AUDIT',
        actual_matlab_execution=True, stage='test', run=str(run), normal_residual_max=threshold,
        artificial_fixture=metadata.get('artificial_fixture') is True,
        inputs=dict(metrics_sha256=input_hashes[str(run/'metrics.csv')],
                    run_metadata_sha256=input_hashes[str(run/'RUN_METADATA.json')],
                    run_status_sha256=input_hashes[str(run/'run_status.json')],
                    experiment_spec_used_sha256=input_hashes[str(run/'experiment_spec_used.json')]),
        scope='Post-hoc numerical reporting only; does not modify scientific acceptance or parameters.',
        counts=counts, solver_summary=solver_summary, selected_route_summary=selected_summary,
        calibration_summary=calibration_summary, issues=issues,
        interpretation_notes=[
            'Original CSV TargetMismatch stores selected Direct mismatch for both routes; raw CSV is unchanged.',
            'Candidate Split mismatch uses split_record.data_residual / sqrt(diagnostics.split_target_squared).',
            'Baseline/Frozen use inherited log.target and legacy route-specific Euclidean residual units; magnitudes are not comparable with Mahalanobis residuals.',
            'Search history includes every searched Direct solve and the final Split solve; selected records are reported separately and are not counted twice in solver_summary.',
            'Initial flags, recorded restart/fallback solvers, and final flags are separate diagnostics. A restart attempt rejected by the solver may not be fully observable from logs.',
            'A Split target mismatch is descriptive: scale selection was frozen to Direct; 0.005 is the existing target-tolerance diagnostic, not a new Split acceptance gate.',
            'Reference/bootstrap optimizer nonpositive exit flags and boundary estimates remain visible; they do not independently prove calibration failure.',
            '56 calibration files and 448 bootstrap estimates, when present, are records with repeated seed namespaces across suites, not 56/448 independent samples.'])
    output.mkdir(parents=True)
    write_csv(output/'case_diagnostics.csv', cases)
    write_csv(output/'solver_records.csv', records)
    write_csv(output/'selected_route_diagnostics.csv', selected)
    write_csv(output/'calibration_estimates.csv', calibration_rows)
    write_csv(output/'calibration_files.csv', calibration_files)
    write_csv(output/'input_manifest.csv', manifest)
    (output/'audit_checks.json').write_text(json.dumps(clean(result), ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    lines = ['# 完整测试求解器与校准后验审计', '',
        f"状态：{result['status']}。仅补充报告，不修改原始CSV、元数据、冻结方法或工程验收。", '',
        f"核对 {len(metrics)} 行指标、{len(cases)} 个方法案例、{len(records)} 条搜索/最终Split求解记录。",
        f"搜索与最终求解记录：最终非零标志 {solver_summary['nonzero_final_flags']}；初始非零标志 {solver_summary['nonzero_initial_flags']}；正规残差超限 {solver_summary['normal_residual_exceedances']}；缺失/非有限记录 {solver_summary['nonfinite_or_missing_records']}；记录的重启或回退 {solver_summary['restart_or_fallback_records']}。", '',
        '原始CSV的TargetMismatch在两条路线均存储Direct值。下表从每条路线自己的残差重算；候选Split使用自己的协方差目标，Baseline/Frozen沿用旧版欧氏单位目标。', '',
        '| 方法 | 路线 | >0.005次数/案例 | 最大实际目标失配 |', '|---|---|---:|---:|']
    for row in selected_summary:
        value = row['max_target_mismatch']
        lines.append(f"|{row['Method']}|{row['Route']}|{row['target_mismatches_gt_0_005']}/{row['cases']}|{value:.6g}|" if value is not None else f"|{row['Method']}|{row['Route']}|不可计算|不可计算|")
    lines += ['', 'Split尺度由冻结的Direct差异原则选择，因此Split失配只是诊断，不新增验收门槛。', '',
        f"校准文件 {len(calibration_files)} 个；参考估计 {calibration_summary['reference_estimates']} 条、bootstrap估计 {calibration_summary['bootstrap_estimates']} 条。模型/校准种子块共 {calibration_summary['unique_model_calibration_seed_blocks']} 个，跨压力组复用；这些计算记录不能当作独立地质或独立校准样本数。", '',
        f"参考估计Q/频率边界：{calibration_summary['reference_q_boundaries']}/{calibration_summary['reference_frequency_boundaries']}；bootstrap：{calibration_summary['bootstrap_q_boundaries']}/{calibration_summary['bootstrap_frequency_boundaries']}。参考/bootstrap优化器非正退出标志：{calibration_summary['reference_optimizer_nonpositive_exitflags']}/{calibration_summary['bootstrap_optimizer_nonpositive_exitflags']}。", '',
        f"缺失/损坏求解日志：{len(missing)}/{corrupt_solver}；缺失/损坏校准日志：{len(missing_calibrations)}/{corrupt_calibration}。所有案例和异常均保留在CSV与audit_checks.json；原始输入字节哈希见input_manifest.csv。", '']
    (output/'NUMERICAL_AUDIT_ZH.md').write_text('\n'.join(lines), encoding='utf-8')
    return result


def self_test():
    class ReportingTests(unittest.TestCase):
        def test_candidate_route_targets_are_distinct(self):
            log = dict(Method='CombinedUncertainty', diagnostics=dict(direct_target_squared=100, split_target_squared=81))
            self.assertEqual(selected_target(log, 'Direct')[0], 10)
            self.assertEqual(selected_target(log, 'Split')[0], 9)
            self.assertAlmostEqual(mismatch(9, selected_target(log, 'Split')[0]), 0)
            self.assertAlmostEqual(mismatch(9, selected_target(log, 'Direct')[0]), .1)

        def test_legacy_uses_inherited_units(self):
            target, units = selected_target(dict(Method='Frozen', target=.03), 'Split')
            self.assertEqual(target, .03)
            self.assertIn('legacy', units)

        def test_initial_failure_is_separate_from_final_success(self):
            d = record_diagnostics(dict(flag=0, initial_flag=1, iterations=40,
                normal_residual=1e-9, data_residual=1, solver='PCG_OriginalPrecision_Restart'), 2e-7)
            self.assertTrue(d['InitialFlagNonzero'])
            self.assertTrue(d['RestartOrFallbackRecorded'])
            self.assertFalse(d['FinalFlagNonzero'])

        def test_null_and_nonfinite_are_retained_as_findings(self):
            d = record_diagnostics(dict(flag=5, initial_flag=0, iterations=20,
                normal_residual=None, data_residual=float('nan')), 2e-7)
            self.assertIn('normal_residual', d['NonfiniteOrMissingFields'])
            self.assertIn('data_residual', d['NonfiniteOrMissingFields'])
            self.assertTrue(d['FinalFlagNonzero'])
            self.assertIsNone(mismatch(None, 10))

        def test_partial_test_is_rejected_without_output(self):
            with tempfile.TemporaryDirectory(prefix='artificial_posthoc_audit_') as temporary:
                root = Path(temporary)
                run = root/'raw'
                run.mkdir()
                (run/'RUN_METADATA.json').write_text(json.dumps(dict(status='RUNNING', stage='test')), encoding='utf-8')
                (run/'run_status.json').write_text(json.dumps(dict(status='RUNNING', stage='test')), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'has not COMPLETED'):
                    run_audit(run, root/'output')
                self.assertFalse((root/'output').exists())

        def test_complete_artificial_case_reports_split_and_all_searches(self):
            # Deliberately artificial schema fixture: no inversion or geological
            # result is generated. The audit must use Split's target, preserve
            # the initial failed search flag, and bind raw-input byte hashes.
            with tempfile.TemporaryDirectory(prefix='artificial_posthoc_audit_') as temporary:
                root = Path(temporary)
                run = root/'raw'
                (run/'solver_logs').mkdir(parents=True)
                (run/'calibration').mkdir()
                meta = dict(status='COMPLETED', stage='test', actual_matlab_execution=True,
                            artificial_fixture=True, rows=2)
                (run/'RUN_METADATA.json').write_text(json.dumps(meta), encoding='utf-8')
                (run/'run_status.json').write_text(json.dumps(meta), encoding='utf-8')
                (run/'experiment_spec_used.json').write_text(json.dumps(dict(algorithm=dict(bootstrap_replicates=8))), encoding='utf-8')
                key = ('ArtificialModel', 'Main', 'SmoothAVO', '101', '201', 'CombinedUncertainty')
                write_csv(run/'metrics.csv', [dict(identity(key), Route=route, Status='SUCCESS', TargetMismatch=.1)
                                              for route in ('Direct', 'Split')])
                record = dict(flag=0, initial_flag=0, iterations=10, normal_residual=1e-9,
                              data_residual=9, solver='PCG_GramWoodbury')
                failed_initial = dict(record, initial_flag=1, solver='PCG_OriginalPrecision_Restart')
                log = dict(identity(key), records=[failed_initial, record, record],
                           direct_record=record, split_record=record, target=10,
                           scale_grid=[.1, 1], scale=1, endpoint=True, residual_mismatch=.1,
                           diagnostics=dict(direct_target_squared=100, split_target_squared=81))
                (run/'solver_logs'/'artificial.json').write_text(json.dumps(log), encoding='utf-8')
                est = dict(effective_Q=40, estimated_frequency=60, minimum_loss=10,
                    optimizer_exitflag=1, optimizer_evaluations=30, q_boundary=False, frequency_boundary=False)
                cal = dict(Model=key[0], Suite=key[1], CalibrationSeed=101, bootstrap_seed=301,
                           estimate=est, bootstrap_estimates=[est]*8)
                (run/'calibration'/'artificial.json').write_text(json.dumps(cal), encoding='utf-8')
                result = run_audit(run, root/'output')
                self.assertEqual(result['solver_summary']['records_total'], 3)
                self.assertEqual(result['solver_summary']['nonzero_initial_flags'], 1)
                split = next(r for r in result['selected_route_summary'] if r['Route'] == 'Split')
                self.assertAlmostEqual(split['max_target_mismatch'], 0)
                direct = next(r for r in result['selected_route_summary'] if r['Route'] == 'Direct')
                self.assertAlmostEqual(direct['max_target_mismatch'], .1)
                self.assertEqual(result['calibration_summary']['bootstrap_estimates'], 8)
                self.assertEqual(result['inputs']['metrics_sha256'], hashlib.sha256((run/'metrics.csv').read_bytes()).hexdigest())
                with (root/'output'/'calibration_estimates.csv').open(encoding='utf-8-sig') as stream:
                    self.assertEqual(len(list(csv.DictReader(stream))), 9)
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReportingTests))
    return 0 if result.wasSuccessful() else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if not args.run or not args.output:
        parser.error('--run and --output are required')
    try:
        result = run_audit(args.run, args.output)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print('AUDIT_REFUSED: ' + str(exc), file=sys.stderr)
        return 2
    print(json.dumps({k: result[k] for k in ('status', 'counts', 'solver_summary')}, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
