from pathlib import Path
import urllib.request,urllib.parse,json,hashlib
from concurrent.futures import ThreadPoolExecutor
R=Path(__file__).resolve().parent
def get(url,path,limit=20_000_000):
    path=R/path;path.parent.mkdir(parents=True,exist_ok=True)
    req=urllib.request.Request(url,headers={'User-Agent':'Research-data-verification/1.0'})
    with urllib.request.urlopen(req,timeout=60) as resp:
        b=resp.read(limit+1)
        if len(b)>limit:raise ValueError('Size cap exceeded')
        info=dict(url=url,final_url=resp.url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),http=dict(resp.headers))
    path.write_bytes(b);return dict(file=str(path.relative_to(R)),**info)
tasks=[]
for name in ['Qp_profile_20-200Hz.npz','README','event_catalog.csv','plot_Q_well_logs.py','read_cat_plot_spectra.py','well_config.npz']:
    tasks.append(('https://gdr.openei.org/files/1823/'+name,'cape/'+name))
for label,filename in [('sumatra','sumatra_metadata.txt'),('whataroa','what_aroa_metadata.txt')]:
    j=json.loads((R/filename).read_text(encoding='utf-8'))
    print(label,j['metadata']['title'],flush=True)
    for f in j['files']:
        print(f['key'],f['size'],f['links'],flush=True)
        if label=='sumatra' and f['size']<5_000_000:tasks.append((f['links']['self'],label+'/'+f['key']))
    (R/(label+'_record_summary.json')).write_text(json.dumps(dict(title=j['metadata']['title'],description=j['metadata'].get('description'),license=j['metadata'].get('license'),files=j['files']),ensure_ascii=False,indent=2),encoding='utf-8')
records=[]
def fetch(task):
    try:return get(*task)
    except Exception as e:return dict(url=task[0],file=task[1],error=str(e))
with ThreadPoolExecutor(max_workers=4) as ex:
    for row in ex.map(fetch,tasks):records.append(row);print(json.dumps(row,ensure_ascii=True),flush=True)
(R/'DOWNLOAD_MANIFEST.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
