"""Seal a completed development config and code for MATLAB held-out execution.

No MATLAB startup/connection. A freeze is a reproducibility guard, not a claim
that the candidate passed. Test receives an explicit expected_hash in its spec.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_manifest(config: dict, here: Path) -> list[dict[str, str]]:
    roots = [Path(config['paths']['source_core']), Path(config['paths']['candidate_matlab']), here]
    roots = [(p if p.is_absolute() else here / p).resolve() for p in roots]
    sources: dict[str, str] = {}
    for root in roots:
        if not root.is_dir():
            raise ValueError(f'Missing source directory: {root}')
        files = list(root.rglob('*.m')) + list(root.rglob('*.py'))
        if (root / 'inputs').is_dir():
            files += list((root / 'inputs').glob('*.json'))
        for path in files:
            if '__pycache__' not in path.parts:
                sources[str(path.resolve())] = sha(path)
    return [{'path': path, 'sha256': sources[path]} for path in sorted(sources)]


def seal(specfile: Path, development: Path, output: Path) -> dict:
    if output.exists():
        raise ValueError('Frozen output must be a new directory')
    specfile = specfile.resolve()
    config = json.loads(specfile.read_text(encoding='utf-8-sig'))
    if 'freeze' in config:
        raise ValueError('Input is already frozen; use the original development spec')
    metadata_path = development / 'RUN_METADATA.json'
    metadata = json.loads(metadata_path.read_text(encoding='utf-8-sig'))
    if metadata.get('status') != 'COMPLETED' or metadata.get('stage') != 'development':
        raise ValueError('Actual completed development run is required before freezing')
    if metadata.get('actual_matlab_execution') is not True:
        raise ValueError('Development was not actually executed in MATLAB')
    if not (development / 'metrics.csv').is_file():
        raise ValueError('Completed development metrics are required')
    if metadata.get('frozen_experiment_spec_sha256') != sha(specfile):
        raise ValueError('Spec changed since development; rerun development before freezing')
    here = Path(__file__).resolve().parent
    sources = source_manifest(config, here)
    prior = {str(Path(x['path']).resolve()): x['sha256'] for x in metadata['source_manifest']}
    current = {x['path']: x['sha256'] for x in sources}
    if prior != current:
        missing = sorted(set(prior) ^ set(current))
        changed = sorted(p for p in set(prior) & set(current) if prior[p] != current[p])
        raise ValueError(f'Sources changed since development; rerun first. Missing: {missing}; changed: {changed}')
    output.mkdir(parents=True)
    snapshot = output / 'config_snapshot.json'
    shutil.copyfile(specfile, snapshot)
    config_sha = sha(snapshot)
    identity = config_sha + '\n' + ''.join(f"{x['path']}={x['sha256']}\n" for x in sources)
    expected = hashlib.sha256(identity.encode('utf-8')).hexdigest()
    config['freeze'] = {
        'expected_hash': expected,
        'config_snapshot': str(snapshot.resolve()),
        'config_sha256': config_sha,
        'sources': sources,
        'development_metadata': str(metadata_path.resolve()),
        'development_metadata_sha256': sha(metadata_path),
        'development_metrics_sha256': sha(development / 'metrics.csv'),
        'development_status': metadata['status'],
        'development_actual_matlab_execution': True,
        'development_failed_rows': metadata.get('failed_rows'),
        'post_test_tuning_allowed': False,
        'hash_is_technical_freeze_not_scientific_acceptance': True,
    }
    (output / 'freeze_identity.txt').write_bytes(identity.encode('utf-8'))
    frozen_spec = output / 'frozen_experiment_spec.json'
    frozen_spec.write_text(json.dumps(config, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    record = {'expected_hash': expected, 'frozen_spec_sha256': sha(frozen_spec),
              'frozen_spec': str(frozen_spec.resolve()), 'post_test_tuning_allowed': False}
    (output / 'FREEZE_RECEIPT.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec', type=Path, required=True)
    parser.add_argument('--development-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(seal(args.spec, args.development_run, args.output), ensure_ascii=False))
    except (ValueError, KeyError, OSError) as exc:
        raise SystemExit(f'FREEZE_REFUSED: {exc}')
