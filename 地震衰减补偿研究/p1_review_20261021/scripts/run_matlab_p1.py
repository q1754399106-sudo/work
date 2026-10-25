"""Connect only to an existing shared ResearchMATLAB session; never start or stop MATLAB."""
from pathlib import Path
import argparse, datetime, io, json, os, sys, time

ROOT = Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--action', choices=['check', 'validate', 'smoke', 'regression', 'replay', 'p1_experiment', 'full'], required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--engine-lib', type=Path)
    a = p.parse_args()
    output = a.output.resolve()
    if output.exists():
        raise SystemExit('Output must be a new directory')
    output.mkdir(parents=True)
    status = {
        'action': a.action,
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'session_required': 'ResearchMATLAB',
        'started_new_instance': False
    }
    before = None; engine = None; out = io.StringIO(); err = io.StringIO(); start = time.monotonic()
    try:
        lib = a.engine_lib or os.environ.get('MATLAB_ENGINE_PYTHON')
        if lib and Path(lib).exists():
            sys.path.insert(0, str(lib))
        try:
            import matlab.engine
        except ModuleNotFoundError as exc:
            status['discovered_sessions'] = []
            status['status'] = 'BLOCKED_ENGINE'
            status['reason'] = f'Official matlab.engine module not present in current runtime ({exc}); no new instance started.'
            return 2
        sessions = matlab.engine.find_matlab()
        status['discovered_sessions'] = list(sessions)
        if 'ResearchMATLAB' not in sessions:
            status['status'] = 'BLOCKED_ENGINE'
            status['reason'] = 'Existing ResearchMATLAB was not discoverable; no new instance started.'
            return 2
        engine = matlab.engine.connect_matlab('ResearchMATLAB')
        def snapshot():
            return {
                'path': engine.path(),
                'pwd': engine.pwd(),
                'rng': engine.feval('jsonencode', engine.rng()),
                'workspace': engine.eval("jsonencode(evalin('base','whos'))"),
                'figures': engine.eval("numel(findall(groot,'Type','figure'))")
            }
        before = snapshot()
        status.update(
            session=engine.feval('matlab.engine.engineName'),
            version=engine.version(),
            sqrt144=engine.sqrt(144.0)
        )
        assert status['session'] == 'ResearchMATLAB' and status['sqrt144'] == 12
        if a.action != 'check':
            target = output / 'matlab'
            engine.addpath(str(ROOT / 'core_p1'), nargout=0)
            engine.rehash(nargout=0)
            assert Path(engine.which('run_p1_research')).resolve() == ROOT / 'core_p1/run_p1_research.m'
            engine.run_p1_research(a.action, str(target), nargout=0, stdout=out, stderr=err)
        status['status'] = 'COMPLETED'
        return 0
    except Exception as exc:
        status['status'] = 'FAILED'
        status['error'] = repr(exc)
        return 1
    finally:
        if before is not None:
            try:
                engine.path(before['path'], nargout=0)
                after = snapshot()
                status['preserved'] = {k: after[k] == v for k, v in before.items()}
                if not all(status['preserved'].values()):
                    status['status'] = 'FAILED_SESSION_PRESERVATION'
            except Exception as exc:
                status['restore_error'] = repr(exc)
                status['status'] = 'FAILED_SESSION_PRESERVATION'
        status['elapsed_seconds'] = time.monotonic() - start
        (output / 'connection.json').write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding='utf-8')
        (output / 'matlab_stdout.log').write_text(out.getvalue(), encoding='utf-8')
        (output / 'matlab_stderr.log').write_text(err.getvalue(), encoding='utf-8')
        print(json.dumps(status, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    raise SystemExit(main())
