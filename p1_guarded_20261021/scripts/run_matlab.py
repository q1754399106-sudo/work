"""Portable hidden MATLAB batch runner; no Engine/SciPy/missing twin dependency."""
from pathlib import Path
import argparse, datetime, hashlib, json, shutil, subprocess, time

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--action',choices=['validate','smoke','guard_trial','full'],default='guard_trial')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--matlab',default=shutil.which('matlab'))
    parser.add_argument('--timeout',type=float,default=1800)
    a=parser.parse_args()
    if not a.matlab: parser.error('MATLAB was not found; supply --matlab PATH')
    output=a.output.resolve()
    if output.exists(): parser.error('Output must be a new directory')
    output.mkdir(parents=True)
    quote=lambda p:p.as_posix().replace("'","''")
    body=f"addpath('{quote(ROOT/'core_p1')}'); "
    if a.action in ('smoke','guard_trial'):
        body+=f"run_p1_research('validate','{quote(output/'validation')}'); "
    body+=f"run_p1_research('{a.action}','{quote(output/'matlab')}');"
    command=[a.matlab,'-wait','-batch',body]
    status={'action':a.action,'command':command,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'status':'RUNNING','timeout_seconds':a.timeout,'source_hashes':{
                str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted((ROOT/'core_p1').rglob('*')) if p.is_file()}}
    start=time.monotonic()
    print(str(output),flush=True)
    try:
        with (output/'stdout.log').open('wb') as out,(output/'stderr.log').open('wb') as err:
            flags=getattr(subprocess,'CREATE_NO_WINDOW',0)
            proc=subprocess.Popen(command,stdout=out,stderr=err,creationflags=flags)
            status['pid']=proc.pid
            (output/'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
            try:
                status['exit_code']=proc.wait(timeout=a.timeout)
                status['status']='COMPLETED' if status['exit_code']==0 else 'FAILED'
            except subprocess.TimeoutExpired:
                print('Hard timeout reached; terminating this batch process tree.',flush=True)
                subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True)
                status['status']='HARD_TIMEOUT'
    finally:
        status['elapsed_seconds']=time.monotonic()-start
        (output/'execution.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
        print(json.dumps({k:v for k,v in status.items() if k!='source_hashes'},indent=2),flush=True)
    return 0 if status['status']=='COMPLETED' else 1

if __name__=='__main__': raise SystemExit(main())
