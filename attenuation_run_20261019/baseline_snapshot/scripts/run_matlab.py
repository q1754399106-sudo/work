"""Connect only to an existing shared session; never start or stop MATLAB."""
from pathlib import Path
import argparse, datetime, io, json, os, sys, time

ROOT = Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--action', choices=['check','validate','regression','replay'], required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--engine-lib', type=Path)
    a = p.parse_args()
    output = a.output.resolve()
    if output.exists(): raise SystemExit('Output must be a new directory')
    output.mkdir(parents=True)
    status = {'action':a.action, 'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'session_required':'ResearchMATLAB','started_new_instance':False}
    before = None; engine = None; out = io.StringIO(); err = io.StringIO(); start = time.monotonic()
    try:
        lib = a.engine_lib or os.environ.get('MATLAB_ENGINE_PYTHON')
        if lib: sys.path.insert(0,str(lib))
        import matlab.engine
        sessions = matlab.engine.find_matlab()
        status['discovered_sessions'] = list(sessions)
        if 'ResearchMATLAB' not in sessions:
            status['status']='BLOCKED_ENGINE'
            status['reason']='Existing ResearchMATLAB was not discoverable; no new instance started.'
            return 2
        engine = matlab.engine.connect_matlab('ResearchMATLAB')
        def snapshot():
            return {'path':engine.path(),'pwd':engine.pwd(),
                    'rng':engine.feval('jsonencode',engine.rng()),
                    'workspace':engine.eval("jsonencode(evalin('base','whos'))"),
                    'figures':engine.eval("numel(findall(groot,'Type','figure'))")}
        before = snapshot()
        status.update(session=engine.feval('matlab.engine.engineName'), version=engine.version(),
                      sqrt144=engine.sqrt(144.0))
        assert status['session']=='ResearchMATLAB' and status['sqrt144']==12
        if a.action!='check':
            # run_research requires a directory that does not yet exist.
            target = output/'matlab'
            engine.addpath(str(ROOT/'core'),nargout=0)
            engine.rehash(nargout=0)
            assert Path(engine.which('run_research')).resolve()==ROOT/'core/run_research.m'
            engine.run_research(a.action,str(target),nargout=0,stdout=out,stderr=err)
        status['status']='COMPLETED'
        return 0
    except Exception as exc:
        status['status']='FAILED'; status['error']=repr(exc)
        return 1
    finally:
        if before is not None:
            try:
                engine.path(before['path'],nargout=0)
                after = snapshot()
                status['preserved']={k:after[k]==v for k,v in before.items()}
                if not all(status['preserved'].values()):
                    status['status']='FAILED_SESSION_PRESERVATION'
            except Exception as exc:
                status['restore_error']=repr(exc); status['status']='FAILED_SESSION_PRESERVATION'
        status['elapsed_seconds']=time.monotonic()-start
        (output/'connection.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
        (output/'matlab_stdout.log').write_text(out.getvalue(),encoding='utf-8')
        (output/'matlab_stderr.log').write_text(err.getvalue(),encoding='utf-8')
        print(json.dumps(status,ensure_ascii=False,indent=2))
        # No quit/exit on the engine: the user owns the shared MATLAB session.

if __name__=='__main__':
    result=main()
    # A restoration failure must never produce a successful exit code.
    args=sys.argv
    if '--output' in args:
        log=Path(args[args.index('--output')+1])/'connection.json'
        if log.exists() and json.loads(log.read_text(encoding='utf-8'))['status'].startswith('FAILED'):
            result=1
    raise SystemExit(result)
