from pathlib import Path
import hashlib, json

ROOT=Path(__file__).resolve().parents[1]
def main():
    records=json.loads((ROOT/'PACKAGE_MANIFEST.json').read_text(encoding='utf-8'))['files']
    errors=[]
    for r in records:
        p=ROOT/r['path']
        if not p.is_file(): errors.append('missing '+r['path'])
        elif hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']: errors.append('changed '+r['path'])
    result={'status':'PASS' if not errors else 'FAIL','files_checked':len(records),'errors':errors}
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return 1 if errors else 0
if __name__=='__main__': raise SystemExit(main())
