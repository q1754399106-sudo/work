"""One-shot existing-session test -> strict audit -> figures; no scheduling."""
from pathlib import Path
import argparse
import datetime as dt
import json
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec',type=Path,required=True)
    parser.add_argument('--tag',default='test_v2')
    args=parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]+',args.tag):
        raise SystemExit('Invalid run tag')
    experiment=ROOT/'experiments'/args.tag
    audit=ROOT/'experiments'/(args.tag+'_strict')
    figures=ROOT/'experiments'/(args.tag+'_figures')
    receipt=ROOT/'experiments'/(args.tag+'_PIPELINE.json')
    if any(p.exists() for p in [experiment,audit,figures,receipt]):
        raise SystemExit('All run, audit, figure and receipt paths must be fresh')
    record=dict(status='RUNNING',started_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                spec=str(args.spec.resolve()),test_run=str(experiment),evaluation=str(audit),
                figures=str(figures),thread_limit_override=None,
                note='Only existing ResearchMATLAB; frozen protocol; all failed rows retained.')
    def save():
        receipt.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    save()
    stages=[('test',[sys.executable,str(ROOT/'candidate/run_candidate.py'),'--action','test',
                     '--output',str(experiment),'--spec',str(args.spec.resolve())]),
            ('strict_audit',[sys.executable,str(ROOT/'candidate/audit/strict_evaluate.py'),
                             '--metrics',str(experiment/'matlab/metrics.csv'),
                             '--spec',str(args.spec.resolve()),
                             '--metadata',str(experiment/'matlab/RUN_METADATA.json'),
                             '--output',str(audit)]),
            ('figures',[sys.executable,str(ROOT/'candidate/visualization/render_results.py'),
                        '--metrics',str(experiment/'matlab/metrics.csv'),
                        '--examples',str(experiment/'matlab/examples'),
                        '--worst-cases',str(experiment/'matlab/worst_cases'),
                        '--output',str(figures)])]
    for stage,command in stages:
        record['active_stage']=stage
        save()
        print('PIPELINE stage '+stage+' START',flush=True)
        log=ROOT/'experiments'/(args.tag+'_'+stage+'_process.log')
        with log.open('w',encoding='utf-8') as stream:
            process=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        record.setdefault('stages',{})[stage]={'returncode':process.returncode,'log':str(log)}
        if process.returncode:
            record.update(status='FAILED_'+stage.upper(),finished_utc=dt.datetime.now(dt.timezone.utc).isoformat())
            save()
            print(json.dumps(record,ensure_ascii=False,indent=2),flush=True)
            return process.returncode
        print('PIPELINE stage '+stage+' COMPLETED',flush=True)
    record.update(status='COMPLETED_REQUIRES_VISUAL_QA_AND_DELIVERY',
                  finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(),active_stage=None)
    record['engineering_result']=json.loads((audit/'acceptance.json').read_text(encoding='utf-8'))['status']
    save()
    print(json.dumps(record,ensure_ascii=False,indent=2),flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
