from pathlib import Path
import json,hashlib,datetime
ROOT=Path(__file__).resolve().parent
paths=[ROOT/'SPEC.json',ROOT/'PLAN_ZH.md']+list((ROOT/'core').rglob('*.m'))+list((ROOT/'core/inputs').glob('*.json'))+list((ROOT/'inputs').glob('*.json'))+list((ROOT/'truth').glob('*.npz'))+[ROOT/p for p in ['analyze.py','run_experiment.py','prepare_data.py','prepare_sensitivity.py','vendor/build_model_original.py']]
out=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 note='Frozen after development and before validation scoring; source snapshots retained. Numerical controls and missing-field fallback added after initial development snapshot; estimator unchanged for that dataset.')
dest=ROOT/'FROZEN_CANDIDATE.json';assert not dest.exists();dest.write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(dict(frozen_utc=out['frozen_utc'],files=len(out['hashes'])),indent=2))
