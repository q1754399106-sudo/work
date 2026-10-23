"""Finish this authorized long-running experiment locally, without modifying its source.

Waits for the runner's terminal status, analyzes every retained record, checks hashes,
and packages reviewable artifacts. Does not silently retry a failed numerical run.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,subprocess,sys,time,zipfile

ROOT=Path(__file__).resolve().parents[1]
def save(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(path)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path);a=parser.parse_args();run=a.run.resolve()
    status=run/'FINALIZATION_STATUS.json'
    save(status,dict(status='WAITING_FOR_NUMERICAL_RUN',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    deadline=time.monotonic()+48*3600
    while not (run/'EXECUTION_SUMMARY.json').exists():
        if time.monotonic()>deadline:
            save(status,dict(status='MONITOR_TIMEOUT',reason='48h elapsed; numerical workers were not terminated'));return 2
        time.sleep(30)
    summary=json.loads((run/'EXECUTION_SUMMARY.json').read_text())
    if summary['status']!='COMPLETED':
        save(status,dict(status='NUMERICAL_RUN_FAILED',summary=summary,automatic_retry=False));return 1
    if summary['source_changed_during_run']:
        save(status,dict(status='SOURCE_CHANGED_REVIEW_REQUIRED',files=summary['source_changed_during_run']));return 1
    save(status,dict(status='ANALYZING',rows=summary['rows']))
    with (run/'analysis_stdout.log').open('wb') as out,(run/'analysis_stderr.log').open('wb') as err:
        result=subprocess.run([sys.executable,str(ROOT/'analysis/analyze_run.py'),str(run)],stdout=out,stderr=err)
    if result.returncode:
        save(status,dict(status='ANALYSIS_FAILED',exit_code=result.returncode,automatic_retry=False));return 1
    coverage=json.loads((run/'analysis/coverage.json').read_text())
    acceptance=json.loads((run/'analysis/acceptance_summary.json').read_text())
    if not coverage['complete']:
        save(status,dict(status='COVERAGE_INCOMPLETE',coverage=coverage));return 1
    files=[p for p in sorted(run.rglob('*')) if p.is_file() and '__pycache__' not in p.parts
           and p.name not in ['FINALIZATION_STATUS.json','PACKAGE_MANIFEST.json'] and not p.name.endswith('.tmp')]
    manifest=dict(files=[dict(path=p.relative_to(run).as_posix(),bytes=p.stat().st_size,
                 sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files])
    manifest['file_count']=len(files);save(run/'PACKAGE_MANIFEST.json',manifest)
    bad=[x['path'] for x in manifest['files'] if hashlib.sha256((run/x['path']).read_bytes()).hexdigest()!=x['sha256']]
    if bad:
        save(status,dict(status='PACKAGE_HASH_FAILURE',files=bad));return 1
    archive=ROOT.parent/(ROOT.name+'_results_'+run.name+'.zip')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[run/'PACKAGE_MANIFEST.json']:z.write(p,p.relative_to(run).as_posix())
        for p in [ROOT/'README_ZH.md',ROOT/'CACHE_EQUIVALENCE_FINAL.json',Path(__file__),ROOT/'analysis/analyze_run.py']:
            z.write(p,'delivery/'+p.name)
    final=dict(status='COMPLETED',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
      experiment_completed=True,coverage=coverage,combined_acceptance=acceptance['CombinedUncertainty'],
      report=str(run/'analysis/RESULTS_ZH.md'),archive=str(archive),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    save(status,final);save(ROOT/'LATEST_RESULT.json',final)
    print(json.dumps(final,ensure_ascii=False,indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
