"""Package code, full tables, verification and frozen-run provenance."""
from pathlib import Path
import shutil,json,hashlib,datetime
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    run=Path(json.loads((ROOT/'LATEST_RUN.json').read_text(encoding='utf-8'))['run'])
    assert json.loads((run/'execution.json').read_text())['status']=='COMPLETED'
    target=ROOT/'deliverable';target.mkdir(exist_ok=False)
    for folder in ['core_p1','scripts','analysis','figures']:
        shutil.copytree(ROOT/folder,target/folder,ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['PLAN_ZH.md','README_ZH.md','RESULTS_ZH.md','LATEST_RUN.json']:
        shutil.copy2(ROOT/name,target/name)
    snapshots=target/'evidence';snapshots.mkdir()
    for rr in sorted((ROOT/'runs').iterdir()):
        if not rr.is_dir():continue
        dst=snapshots/rr.name;dst.mkdir()
        if (rr/'source').exists():shutil.copytree(rr/'source',dst/'source')
        for name in ['execution.json','SOURCE_MANIFEST.json','stdout.log','stderr.log']:
            if (rr/name).exists():shutil.copy2(rr/name,dst/name)
        if (rr/'matlab').exists():
            (dst/'matlab').mkdir()
            for p in (rr/'matlab').iterdir():
                if p.suffix in ('.json','.csv'):shutil.copy2(p,dst/'matlab'/p.name)
    current_numeric=['joint_select.m','joint_physical_ensemble.m','p1_solve_branch.m']
    for name in current_numeric:assert sha(ROOT/'core_p1/matlab_src'/name)==sha(run/'source/core_p1/matlab_src'/name)
    # Compare copied dependencies with the untouched prior project.
    baseline=ROOT.parent/'p1_improved_20261022/core_p1'
    original_differences=[]
    for p in baseline.rglob('*'):
        if p.is_file():
            rel=p.relative_to(baseline);q=ROOT/'core_p1'/rel
            if not q.exists() or sha(p)!=sha(q):original_differences.append(rel.as_posix())
    assert original_differences==['matlab_src/p1_solve_branch.m'],original_differences
    manifest={p.relative_to(target).as_posix():sha(p) for p in sorted(target.rglob('*')) if p.is_file()}
    (target/'MANIFEST.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'prior_dependency_difference':original_differences,'files':manifest},ensure_ascii=False,indent=2),encoding='utf-8')
    archive=shutil.make_archive(str(ROOT.parent/'joint_refactor_20261022_results'),'zip',target)
    print(json.dumps({'archive':archive,'bytes':Path(archive).stat().st_size,'sha256':sha(Path(archive)),'files':len(manifest)},indent=2))
if __name__=='__main__':main()
