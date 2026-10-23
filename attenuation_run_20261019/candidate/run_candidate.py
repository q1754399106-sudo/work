"""Run candidate functions through the user's existing shared MATLAB only."""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import io
import json
import sys
import time

ROOT = Path(__file__).resolve().parent
RUN = ROOT.parent
ENGINE_LIB = Path('C:/Users/qq175/Documents/Codex/2026-10-17/new-chat/matlab_engine_runtime/lib')


class LiveLog(io.StringIO):
    def __init__(self, path):
        super().__init__()
        self.handle = path.open('w', encoding='utf-8', buffering=1)
    def write(self, text):
        super().write(text)
        self.handle.write(text)
        self.handle.flush()
        print(text, end='', flush=True)
        return len(text)
    def flush(self):
        self.handle.flush()
    def close(self):
        if not self.handle.closed:
            self.handle.close()
        super().close()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--action', required=True, choices=['selftest', 'development', 'test', 'checkcode', 'equivalence'])
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--spec', type=Path, default=ROOT/'experiment'/'experiment_spec.json')
    parser.add_argument('--engine-lib', type=Path, default=ENGINE_LIB)
    parser.add_argument('--threads', type=int, help='Temporary computational-thread limit, restored on return')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    status = dict(action=args.action, started_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                  existing_session_required='ResearchMATLAB', started_new_instance=False,
                  original_package_modified=False)
    started = time.monotonic()
    stdout = LiveLog(output/'matlab_stdout.log')
    stderr = LiveLog(output/'matlab_stderr.log')
    engine = None
    before = None
    result = 1
    try:
        sys.path.insert(0, str(args.engine_lib))
        import matlab.engine
        sessions = matlab.engine.find_matlab()
        status['discovered_sessions'] = list(sessions)
        if 'ResearchMATLAB' not in sessions:
            status.update(status='BLOCKED_ENGINE', reason='No new MATLAB instance was started.')
            return 2
        engine = matlab.engine.connect_matlab('ResearchMATLAB')
        def snapshot():
            return dict(path=engine.path(), pwd=engine.pwd(), rng=engine.feval('jsonencode',engine.rng()),
                        computational_threads=engine.maxNumCompThreads(),
                        workspace=engine.eval("jsonencode(evalin('base','whos'))"),
                        figures=engine.eval("numel(findall(groot,'Type','figure'))"))
        before = snapshot()
        rng_before = engine.rng()
        if args.threads is not None:
            assert args.threads >= 1
            engine.maxNumCompThreads(float(args.threads), nargout=0)
        status['computational_threads_used']=engine.maxNumCompThreads()
        status.update(session=engine.feval('matlab.engine.engineName'), version=engine.version(),
                      sqrt144=engine.sqrt(144.0))
        assert status['session']=='ResearchMATLAB' and status['sqrt144']==12
        folders=[RUN/'baseline_snapshot'/'core'/'matlab_src', ROOT/'matlab', ROOT/'experiment']
        if args.action=='equivalence':
            folders.append(ROOT/'numerical_equivalence')
        for folder in folders:
            engine.addpath(str(folder), nargout=0)
        engine.rehash(nargout=0)
        targets={'p1_compensate':ROOT/'matlab'/'p1_compensate.m',
                 'p1_run_experiment':ROOT/'experiment'/'p1_run_experiment.m',
                 'prepare_models_v2':RUN/'baseline_snapshot'/'core'/'matlab_src'/'prepare_models_v2.m'}
        for name, expected in targets.items():
            actual = engine.which(name)
            assert Path(actual).resolve()==expected.resolve(), f'Unexpected path shadowing: {name}: {actual}'
        status['resolved_sources']={name:str(path) for name,path in targets.items()}
        status['source_hashes']={str(p.relative_to(RUN)):sha(p) for p in ROOT.rglob('*')
                                 if p.is_file() and p.suffix in {'.m','.py','.json'} and '__pycache__' not in p.parts}
        if args.action=='selftest':
            engine.p1_selftest(str(output/'matlab'), nargout=0, stdout=stdout, stderr=stderr)
        elif args.action=='equivalence':
            engine.p1_compare_numerics(str(output/'matlab'), nargout=0, stdout=stdout, stderr=stderr)
        elif args.action=='checkcode':
            diagnostics={}
            for p in ROOT.rglob('*.m'):
                diagnostics[str(p.relative_to(ROOT))] = json.loads(engine.p1_code_diagnostics(str(p)))
            (output/'checkcode.json').write_text(json.dumps(diagnostics,indent=2),encoding='utf-8')
        else:
            status['spec_sha256']=sha(args.spec)
            engine.p1_run_experiment(args.action, str(output/'matlab'), str(args.spec.resolve()),
                                     nargout=0, stdout=stdout, stderr=stderr)
        status['status']='COMPLETED'
        result = 0
    except Exception as exc:
        status.update(status='FAILED', error=repr(exc))
        result = 1
    finally:
        if before is not None:
            try:
                engine.path(before['path'], nargout=0)
                engine.cd(before['pwd'], nargout=0)
                engine.rng(rng_before, nargout=0)
                if args.threads is not None:
                    engine.maxNumCompThreads(before['computational_threads'], nargout=0)
                after = snapshot()
                status['preserved']={key:after[key]==value for key,value in before.items()}
                status['preservation_scope']='Path, pwd, RNG and computational-thread count exact; base variable metadata and figure count checked; candidate functions never write base or figures.'
                if not all(status['preserved'].values()):
                    status['status']='FAILED_SESSION_PRESERVATION'
                    result=1
            except Exception as exc:
                status.update(status='FAILED_SESSION_PRESERVATION',restore_error=repr(exc))
                result=1
        status['elapsed_seconds']=time.monotonic()-started
        (output/'connection.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
        stdout.close()
        stderr.close()
        concise={k:v for k,v in status.items() if k not in {'source_hashes','resolved_sources'}}
        print(json.dumps(concise,ensure_ascii=False,indent=2))
    return result


if __name__=='__main__':
    raise SystemExit(main())
