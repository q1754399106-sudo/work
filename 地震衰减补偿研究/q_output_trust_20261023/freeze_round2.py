from pathlib import Path
import json,datetime,hashlib
ROOT=Path(__file__).resolve().parent
assert (ROOT/'SELECTED.json').exists() and (ROOT/'inputs/validation_generated.json').exists()
paths=list((ROOT/'core').rglob('*.m'))+list((ROOT/'core/inputs').glob('*.json'))+list((ROOT/'inputs').glob('*.json'))+list((ROOT/'truth').glob('*.npz'))+[ROOT/p for p in ['SPEC.json','PLAN_ZH.md','SELECTED.json','run_round2.py','analyze_round2.py','prepare_data.py','prepare_sensitivity.py','summarize_sensitivity.py','vendor/build_model_original.py']]
out=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'FROZEN_CANDIDATE.json';assert not dest.exists();dest.write_text(json.dumps(out,indent=2),encoding='utf-8');print('Frozen',len(paths),'files',out['frozen_utc'])
