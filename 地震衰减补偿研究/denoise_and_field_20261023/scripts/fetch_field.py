from pathlib import Path
import urllib.request,json,hashlib,time,concurrent.futures
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'field_data';OUT.mkdir(exist_ok=True)
sources=[dict(name='wz.25.H',url='https://ahay.org/data/wz/wz.25.H',reference='https://ahay.org/RSF/book/tccs/orthogroll/field.html',description='Yilmaz shot record used in published field ground-roll experiments; inspect header and upstream provenance before using.'),dict(name='L23535.SGY',url='https://ahay.org/data/alaska/L23535.SGY',reference='https://ahay.org/RSF/book/data/alaska/line31-81-m8r.html',description='Madagascar mirror of USGS Alaska NPRA line 31-81 raw seismic; one original SEG-Y file.')]
def fetch(s):
 start=time.monotonic();p=OUT/s['name'];r=dict(s)
 try:
  request=urllib.request.Request(s['url'],headers={'User-Agent':'Research-data-validation/1.0'})
  with urllib.request.urlopen(request,timeout=45) as response,p.open('wb') as f:
   r['headers']=dict(response.headers);r['status']=response.status
   n=0
   while True:
    b=response.read(1024*1024)
    if not b:break
    n+=len(b)
    if n>150_000_000:raise RuntimeError('Safety size cap 150MB exceeded; partial retained')
    f.write(b)
  r.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-start,downloaded=True)
 except Exception as e:r.update(downloaded=False,error=str(e),elapsed_seconds=time.monotonic()-start)
 (OUT/(s['name']+'.download.json')).write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(r,ensure_ascii=True),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(2) as executor:results=list(executor.map(fetch,sources))
(OUT/'download_manifest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
