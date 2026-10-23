from pathlib import Path
import collections, csv, hashlib, json, statistics

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'review_run'
manifest = json.loads((ROOT / 'PACKAGE_MANIFEST.json').read_text(encoding='utf-8'))
checks = []
for f in manifest['files']:
    p = ROOT / f['path']
    checks.append({'path': f['path'], 'exists': p.is_file(),
                   'hash_matches': p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() == f['sha256'],
                   'size_matches': p.is_file() and p.stat().st_size == f['bytes']})
result = {'manifest_declared_count': manifest['file_count'], 'entries_checked': len(checks),
          'mismatches': [c for c in checks if not all(c[k] for k in ('exists','hash_matches','size_matches'))]}
with (ROOT / 'results/metrics_p1_all_models_multical.csv').open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
result['historical_file_rows'] = len(rows)
result['split_branch_summary'] = {}
for method in sorted({r['Method'] for r in rows}):
    rr = [r for r in rows if r['Method'] == method and r['Route'] == 'Split']
    result['split_branch_summary'][method] = {'observations': len(rr),
        'failures_error_gt_1': sum(float(r['RelativeError']) > 1 for r in rr)}
result['calibration_observation_counts'] = []
for model in sorted({r['Model'] for r in rows}):
    rr = [r for r in rows if r['Model']==model and r['Method']=='CombinedUncertainty' and r['Route']=='Split']
    result['calibration_observation_counts'].append({'model':model,
        'counts':dict(collections.Counter((r['Suite']+'/'+r['CalibrationSeed']) for r in rr))})
for name in ('historical_recheck', 'all_models_recheck'):
    r = json.loads((OUT / name / 'acceptance.json').read_text(encoding='utf-8'))
    result[name] = {k:r[k] for k in ('status','pairs','equal_model_main_improvement')}
    result[name]['failed_checks'] = [c for c in r['checks'] if not c['pass_']]
(OUT / 'independent_audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=False, indent=2))
