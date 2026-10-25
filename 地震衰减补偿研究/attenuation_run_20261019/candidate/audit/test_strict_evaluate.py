"""Artificial fixtures exercise audit logic; these are NOT scientific results."""
from pathlib import Path
import copy
import csv
import json
import tempfile
import unittest

import strict_evaluate as audit


class StrictAuditFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture_spec = dict(
            status='ARTIFICIAL_UNIT_TEST_NOT_AN_EXPERIMENT',
            development_model_ids=['ArtificialDevA', 'ArtificialDevB'],
            test_model_ids=['ArtificialTest' + str(x) for x in range(4)],
            candidate_branches=['NoiseUncertaintyOnly', 'QModelErrorOnly', 'CombinedUncertainty'],
            controls=[], scenario_names=['SmoothAVO', 'FaultSection', 'WeakSection', 'PinchoutSection'],
            snr_db=[18, 12], suite_snr_db={'Main': [18], 'LowerSNR': [12]},
            required_suites=['Main', 'LowerSNR'], bootstrap_replicates=100, bootstrap_seed=19)
        cls.fixture_rows = []
        for model in cls.fixture_spec['test_model_ids']:
            for suite, snr in [('Main', 18), ('LowerSNR', 12)]:
                for calibration in range(3):
                    for seed in range(6):
                        for scenario in cls.fixture_spec['scenario_names']:
                            for method in ['Baseline', 'Frozen'] + cls.fixture_spec['candidate_branches']:
                                for route in ['Direct', 'Split']:
                                    scale = .8 if method in cls.fixture_spec['candidate_branches'] else 1.
                                    row = dict(Model=model, Suite=suite, SNRdB=str(snr),
                                        CalibrationSeed=str(100 + calibration), Seed=str(200 + seed),
                                        Scenario=scenario, Method=method, Route=route,
                                        PriorHash='a' * 64, SharedInformationHash='b' * 64,
                                        NormalResidual='1e-10', FinalFlag='0', Status='COMPLETED',
                                        FailureMessage='', ScaleEndpoint='0', TargetMismatch='0')
                                    row.update({metric: str((.8 if metric == 'RelativeError' else .1) * scale)
                                                if metric in audit.applicable(scenario) else '' for metric in audit.METRICS})
                                    cls.fixture_rows.append(row)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='artificial_audit_fixture_')
        self.root = Path(self.temporary.name)
        self.source = self.root / 'artificial_source.py'
        self.source.write_text('# Artificial test input only. No inversion implemented.\n', encoding='utf-8')
        self.spec = self.root / 'frozen_spec.json'
        self.spec.write_text(json.dumps(self.fixture_spec, indent=2), encoding='utf-8')
        self.metadata = self.root / 'RUN_METADATA.json'
        self.metadata_record = dict(artificial_fixture=True, actual_matlab_execution=False,
            status='ARTIFICIAL_TEST', all_failures_retained=True,
            acceptance_thresholds_frozen_before_test=True,
            frozen_experiment_spec_sha256=audit.sha256(self.spec),
            source_manifest=[dict(path=str(self.source), sha256=audit.sha256(self.source))],
            truth_isolation_audit={'scope': 'Artificial fixture; no geological data or truth'},
            new_geological_holdout={'status': 'NOT_APPLICABLE_ARTIFICIAL_FIXTURE'})
        self.metrics = self.root / 'artificial_metrics.csv'

    def tearDown(self):
        self.temporary.cleanup()

    def run_fixture(self, rows=None, metadata=None, spec=None):
        rows = self.fixture_rows if rows is None else rows
        if spec is not None:
            self.spec.write_text(json.dumps(spec, indent=2), encoding='utf-8')
            self.metadata_record['frozen_experiment_spec_sha256'] = audit.sha256(self.spec)
        record = self.metadata_record if metadata is None else metadata
        self.metadata.write_text(json.dumps(record, indent=2), encoding='utf-8')
        if rows is not False:
            with self.metrics.open('w', encoding='utf-8', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(self.fixture_rows[0]))
                writer.writeheader()
                writer.writerows(rows)
        return audit.evaluate(self.metrics, self.spec, self.metadata, self.root / 'result')

    def test_complete_artificial_fixture_passes_without_becoming_science(self):
        result = self.run_fixture()
        self.assertEqual(result['status'], 'ARTIFICIAL_FIXTURE_PASS')
        self.assertFalse(self.metadata_record['actual_matlab_execution'])
        self.assertEqual(result['observations'], 576)
        branch = result['branches']['CombinedUncertainty']['routes']['Split']
        self.assertAlmostEqual(branch['equal_model_main_improvement'], .2)
        interval = branch['hierarchical_uncertainty']['relative_improvement_interval']
        self.assertAlmostEqual(interval[0], .2)
        self.assertAlmostEqual(interval[1], .2)
        self.assertTrue((self.root / 'result' / 'paired_differences.csv').is_file())

    def test_no_csv_reports_blocked_engine_without_fabrication(self):
        self.metadata_record.update(status='BLOCKED_ENGINE', artificial_fixture=False)
        result = self.run_fixture(rows=False)
        self.assertEqual(result['status'], 'BLOCKED_ENGINE')
        self.assertNotIn('branches', result)

    def test_empty_csv_reports_not_run(self):
        self.assertEqual(self.run_fixture(rows=[])['status'], 'NOT_RUN')

    def test_duplicate_observation_method_route_rejected(self):
        result = self.run_fixture(self.fixture_rows + [dict(self.fixture_rows[0])])
        self.assertEqual(result['status'], 'INVALID_INPUT')
        self.assertIn('Duplicate case', result['reason'])

    def test_missing_calibration_layer_rejected(self):
        rows = [r for r in self.fixture_rows if r['CalibrationSeed'] != '102']
        self.assertEqual(self.run_fixture(rows)['status'], 'INVALID_INPUT')

    def test_missing_observation_layer_rejected(self):
        rows = [r for r in self.fixture_rows if r['Seed'] != '205']
        self.assertEqual(self.run_fixture(rows)['status'], 'INVALID_INPUT')

    def test_local_coverage_required_for_each_model_snr(self):
        rows = [r for r in self.fixture_rows if not (r['Model'] == 'ArtificialTest0'
            and r['SNRdB'] == '12' and r['Scenario'] == 'WeakSection')]
        result = self.run_fixture(rows)
        self.assertEqual(result['status'], 'INVALID_INPUT')
        self.assertIn('Missing model/SNR', result['reason'])

    def test_matched_prior_and_information_hashes_enforced(self):
        rows = copy.deepcopy(self.fixture_rows)
        rows[0]['PriorHash'] = 'c' * 64
        result = self.run_fixture(rows)
        self.assertEqual(result['status'], 'INVALID_INPUT')
        self.assertIn('prior hash differs', result['reason'])

    def test_shared_information_hash_across_all_methods(self):
        rows = copy.deepcopy(self.fixture_rows)
        rows[0]['SharedInformationHash'] = 'c' * 64
        result = self.run_fixture(rows)
        self.assertEqual(result['status'], 'INVALID_INPUT')
        self.assertIn('Information differs', result['reason'])

    def test_main_gain_cannot_hide_lower_snr_harm(self):
        rows = copy.deepcopy(self.fixture_rows)
        for row in rows:
            if row['Method'] == 'CombinedUncertainty' and row['SNRdB'] == '12':
                row['RelativeError'] = '.82'
        result = self.run_fixture(rows)
        branch = result['branches']['CombinedUncertainty']
        self.assertEqual(branch['status'], 'FAIL')
        self.assertGreater(branch['routes']['Split']['equal_model_main_improvement'], .05)
        self.assertTrue(any(not x['passed'] and 'LowerSNR' in x['name'] and 'SNR12' in x['name']
                            for x in branch['routes']['Split']['checks']))

    def test_global_gain_cannot_hide_local_weak_harm(self):
        rows = copy.deepcopy(self.fixture_rows)
        for row in rows:
            if row['Method'] == 'CombinedUncertainty' and row['Scenario'] == 'WeakSection':
                row['WeakEventError'] = '.2'
        result = self.run_fixture(rows)
        self.assertEqual(result['branches']['CombinedUncertainty']['status'], 'FAIL')

    def test_failed_direct_nan_row_retained_and_branch_fails(self):
        rows = copy.deepcopy(self.fixture_rows)
        row = next(r for r in rows if r['Method'] == 'CombinedUncertainty' and r['Route'] == 'Direct')
        row.update(Status='FAILED', FailureMessage='ARTIFICIAL_SOLVER_FAILURE', FinalFlag='1', NormalResidual='NaN')
        row.update({metric: 'NaN' for metric in audit.METRICS})
        result = self.run_fixture(rows)
        self.assertEqual(result['branches']['CombinedUncertainty']['status'], 'FAIL')
        self.assertEqual(len(result['failed_rows']), 1)
        self.assertEqual(result['rows'], len(rows))

    def test_successful_nonfinite_metric_rejected(self):
        rows = copy.deepcopy(self.fixture_rows)
        rows[0]['RelativeError'] = 'NaN'
        self.assertEqual(self.run_fixture(rows)['status'], 'INVALID_INPUT')

    def test_preparation_failure_retains_both_routes_with_unavailable_prior(self):
        rows = copy.deepcopy(self.fixture_rows)
        observation = tuple(rows[0][key] for key in audit.KEYS)
        for row in rows:
            if row['Method'] == 'CombinedUncertainty' and tuple(row[key] for key in audit.KEYS) == observation:
                row.update(Status='FAILED', FailureMessage='ARTIFICIAL_CALIBRATION_PREPARATION_FAILURE',
                           PriorHash='FAILED', SharedInformationHash='ARTIFICIAL_INPUT_HASH_ONLY',
                           FinalFlag='1', NormalResidual='NaN')
                row.update({metric: 'NaN' for metric in audit.METRICS})
        result = self.run_fixture(rows)
        self.assertEqual(result['branches']['CombinedUncertainty']['status'], 'FAIL')
        self.assertEqual(len(result['failed_rows']), 2)

    def test_frozen_spec_hash_required(self):
        self.metadata_record['frozen_experiment_spec_sha256'] = '0' * 64
        result = self.run_fixture()
        self.assertEqual(result['status'], 'INVALID_INPUT')
        self.assertIn('specification hash', result['reason'])

    def test_current_sources_must_match_manifest(self):
        self.source.write_text('# Changed artificial source\n', encoding='utf-8')
        result = self.run_fixture()
        self.assertEqual(result['status'], 'INVALID_INPUT')
        self.assertIn('Source hash mismatch', result['reason'])

    def test_scientific_metadata_cannot_substitute_a_boolean_for_freeze_envelope(self):
        record = dict(self.metadata_record, artificial_fixture=False,
                      actual_matlab_execution=True, status='COMPLETED')
        result = self.run_fixture(metadata=record)
        self.assertEqual(result['status'], 'INVALID_INPUT')
        self.assertIn('frozen source/config envelope', result['reason'])

    def test_source_configuration_freeze_identity_checked(self):
        # Only audit primitives are tested here. No scientific metrics or pass
        # artifact is emitted from these deliberately artificial metadata.
        original = copy.deepcopy(self.fixture_spec)
        snapshot = self.root / 'ARTIFICIAL_CONFIG_SNAPSHOT.json'
        snapshot.write_text(json.dumps(original, indent=2), encoding='utf-8')
        sources = self.metadata_record['source_manifest']
        freeze = dict(config_snapshot=str(snapshot), config_sha256=audit.sha256(snapshot),
            sources=sources, development_status='COMPLETED', development_actual_matlab_execution=True)
        identity = freeze['config_sha256'] + '\n' + ''.join(
            r['path'] + '=' + r['sha256'] + '\n' for r in sources)
        freeze['expected_hash'] = audit.hashlib.sha256(identity.encode('utf-8')).hexdigest()
        frozen = dict(original, freeze=freeze)
        self.spec.write_text(json.dumps(frozen, indent=2), encoding='utf-8')
        record = dict(self.metadata_record, artificial_fixture=False,
                      actual_matlab_execution=True, status='COMPLETED',
                      frozen_experiment_spec_sha256=audit.sha256(self.spec))
        self.assertEqual(audit.check_metadata(record, self.metadata, self.spec, frozen, False), 1)
        freeze['expected_hash'] = '0' * 64
        self.spec.write_text(json.dumps(frozen, indent=2), encoding='utf-8')
        record['frozen_experiment_spec_sha256'] = audit.sha256(self.spec)
        with self.assertRaisesRegex(audit.InvalidInput, 'identity hash mismatch'):
            audit.check_metadata(record, self.metadata, self.spec, frozen, False)

    def test_frozen_stress_subset_accepted_without_claiming_calibration_variance(self):
        spec = copy.deepcopy(self.fixture_spec)
        spec['suite_snr_db']['StressImpulse'] = [18]
        spec['required_suites'].append('StressImpulse')
        spec['suites'] = [dict(name=suite, snr_db=snr, scenario_names=spec['scenario_names'],
                              calibration_count=3, observation_count=6)
                          for suite, snr in [('Main', 18), ('LowerSNR', 12)]]
        spec['suites'].append(dict(name='StressImpulse', snr_db=18,
            scenario_names=['WeakSection'], calibration_count=1, observation_count=6))
        extra = []
        for row in self.fixture_rows:
            if row['Suite'] == 'Main' and row['CalibrationSeed'] == '100' and row['Scenario'] == 'WeakSection':
                changed = dict(row)
                changed['Suite'] = 'StressImpulse'
                extra.append(changed)
        result = self.run_fixture(self.fixture_rows + extra, spec=spec)
        self.assertEqual(result['status'], 'ARTIFICIAL_FIXTURE_PASS')

    def test_four_white_scenes_and_coherent_only_in_frozen_stress_suite(self):
        spec = copy.deepcopy(self.fixture_spec)
        white_scenes = list(spec['scenario_names'])
        spec['scenario_names'].append('CoherentNoiseAVO')
        spec['suite_snr_db']['StressCoherent'] = [18]
        spec['required_suites'].append('StressCoherent')
        spec['suites'] = [dict(name=suite, snr_db=snr, scenario_names=white_scenes,
                              calibration_count=3, observation_count=6)
                          for suite, snr in [('Main', 18), ('LowerSNR', 12)]]
        spec['suites'].append(dict(name='StressCoherent', snr_db=18,
            scenario_names=['CoherentNoiseAVO'], calibration_count=1, observation_count=6))
        extra = []
        for row in self.fixture_rows:
            if row['Suite'] == 'Main' and row['CalibrationSeed'] == '100' and row['Scenario'] == 'SmoothAVO':
                changed = dict(row)
                changed.update(Suite='StressCoherent', Scenario='CoherentNoiseAVO')
                extra.append(changed)
        result = self.run_fixture(self.fixture_rows + extra, spec=spec)
        self.assertEqual(result['status'], 'ARTIFICIAL_FIXTURE_PASS')
        self.assertEqual(result['observations'], 576 + 24)


if __name__ == '__main__':
    unittest.main(verbosity=2)
