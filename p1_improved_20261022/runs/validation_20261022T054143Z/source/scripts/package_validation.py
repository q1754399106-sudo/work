"""Package a terminal, analyzed experiment with per-file hashes and ZIP checks.

Packaging a failed method is permitted and preserves its failures. This command
does not promote the method or claim that the user's research goal is complete.
"""
from pathlib import Path
import argparse,hashlib,json,zipfile
EXCLUDED={'PACKAGE_MANIFEST.json','PACKAGE_VERIFICATION.json'}
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def package(run,output):
 run=run.resolve();output=output.resolve()
 assert not output.exists(),'Use a new archive path; retained results are not overwritten'
 assert not output.is_relative_to(run),'Archive must be outside experiment directory'
 execution=read(run/'EXECUTION_SUMMARY.json');assert execution['status']=='COMPLETED','Numerical tasks must be terminal and completed'
 assert all(x['status']=='COMPLETED' for x in execution['tasks'])
 acceptance=read(run/'analysis/acceptance.json');assert acceptance['coverage']['complete'],'Cannot deliver incomplete coverage'
 assert (run/'analysis/RESULTS_ZH.md').is_file()
 for item in ['RUN_METADATA.json','FROZEN_EXECUTION_PLAN.json','ORIGINAL_REPLAY_AUDIT.json','SOURCE_MANIFEST.json','analysis/crossed_bootstrap.json','analysis/integrity.json']:
  assert (run/item).is_file(),f'Missing deliverable: {item}'
 files=[p for p in sorted(run.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name not in EXCLUDED]
 manifest={p.relative_to(run).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files}
 assert not any(p.suffix=='.tmp' for p in files),'Incomplete driver write detected'
 save(run/'PACKAGE_MANIFEST.json',dict(status=acceptance['status'],entries=manifest,
  scope='All experiment files except package manifest, post-package verification and Python bytecode; source, data, logs, failures and analyses included'))
 with zipfile.ZipFile(output,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files+[run/'PACKAGE_MANIFEST.json']:z.write(p,p.relative_to(run).as_posix())
 with zipfile.ZipFile(output) as z:
  assert z.testzip() is None,'ZIP CRC failure'
  assert set(z.namelist())==set(manifest)|{'PACKAGE_MANIFEST.json'},'ZIP member mismatch'
  for name,item in manifest.items():
   assert hashlib.sha256(z.read(name)).hexdigest()==item['sha256'],f'Archived hash mismatch: {name}'
 # Recheck originals: packaging must not silently capture concurrently changing data.
 changed=[name for name,item in manifest.items() if sha(run/name)!=item['sha256']]
 assert not changed,f'Files changed during packaging: {changed}'
 verification=dict(verified=True,archive=str(output),archive_sha256=sha(output),archive_bytes=output.stat().st_size,
  manifest_files=len(manifest),zip_crc_passed=True,archived_file_hashes_passed=True,original_files_unchanged=True,method_acceptance=acceptance['status'],goal_completion_not_implied=True)
 save(run/'PACKAGE_VERIFICATION.json',verification);return verification
def main():
 p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 print(json.dumps(package(a.run,a.output),ensure_ascii=True,indent=2))
if __name__=='__main__':main()
