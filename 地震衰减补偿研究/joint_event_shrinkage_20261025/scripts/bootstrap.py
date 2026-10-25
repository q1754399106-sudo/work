"""Pin input bytes and vendor the existing reference into this task only."""
from pathlib import Path
import subprocess,hashlib,json
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parent
BASE='ff95a69e31b9eee462d1a8843960a00b829cb47f'
def gitbytes(path,rev=BASE):
    raw=subprocess.run(['git','-C',str(REPO),'show',f'{rev}:{path}'],capture_output=True,check=True).stdout
    if raw.startswith(b'version https://git-lfs.github.com/spec/v1'):
        fields=dict(x.split(' ',1) for x in raw.decode().strip().splitlines());oid=fields['oid'].split(':')[1]
        raw=(REPO/'.git/lfs/objects'/oid[:2]/oid[2:4]/oid).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==oid and len(raw)==int(fields['size'])
    return raw
def main():
    sources=['novelty_evidence_20261025/scripts/common.py','novelty_evidence_20261025/scripts/benchmark.py','novelty_evidence_20261025/scripts/regional_diagnostic.py',
 'round3_20261024/run_combo.m','round3_20261024/source/reflectivity_refit_prepared.m','round3_20261024/source/reflectivity_soft_prepared.m','round3_20261024/source/refit_dictionary.m',
 'round3_20261024/core/coherence_scores.m','round3_20261024/core/score_gather.m','round3_20261024/core/score_weak.m','qscan_predictor_20261025/scripts/full_refit.m']
    manifest={}
    for src in sources:
        raw=gitbytes(src);dst=ROOT/'inputs/source_snapshots'/src;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(raw)
        manifest[src]={'commit':BASE,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
    common=gitbytes(sources[0]).decode('utf-8')
    common=common.replace("str(ROOT/'cache'/'deps')","str(REPO/'novelty_evidence_20261025'/'cache'/'deps')")
    common=common.replace("raw=REPO/f'round3_20261024/cache/{name}_Q+0_{noise}_nominal.mat'", "raw=ROOT/f'inputs/prepared/{name}_Q+0_{noise}_nominal.mat'\n    if not raw.exists():\n        raw.parent.mkdir(parents=True,exist_ok=True)\n        raw.write_bytes(prepared_bytes(name,noise))")
    common += '\n\ndef prepared_bytes(name,noise):\n    return pinned_bytes(f"novelty_evidence_20261025/inputs/prepared/{name}_Q+0_{noise}_nominal.mat")\n\ndef pinned_bytes(path):\n    from bootstrap import gitbytes as read_pinned\n    return read_pinned(path)\n'
    (ROOT/'scripts/common.py').write_text(common,encoding='utf-8')
    (ROOT/'scripts/reference.py').write_bytes(gitbytes(sources[1]))
    (ROOT/'inputs/source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Vendored reference and source snapshots:',len(manifest))
if __name__=='__main__':main()
