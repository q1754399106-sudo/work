"""Read public ZIP central directories and a few members using HTTP ranges."""
from pathlib import Path
import io,json,urllib.request,zipfile,hashlib,csv
R=Path(__file__).resolve().parent
class Remote(io.RawIOBase):
    def __init__(self,url):
        self.url=url;self.pos=0;self.bytes=0;self.requests=[]
        req=urllib.request.Request(url,method='HEAD',headers={'User-Agent':'Research-data-verification/1.0'})
        with urllib.request.urlopen(req,timeout=45) as resp:self.size=int(resp.headers['Content-Length'])
    def seekable(self):return True
    def readable(self):return True
    def tell(self):return self.pos
    def seek(self,off,whence=0):
        self.pos=off if whence==0 else self.pos+off if whence==1 else self.size+off
        if self.pos<0:raise ValueError('negative seek')
        return self.pos
    def read(self,n=-1):
        if n<0:n=self.size-self.pos
        n=min(n,self.size-self.pos)
        if n<=0:return b''
        if n>12_000_000 or self.bytes+n>35_000_000:raise ValueError('Remote transfer cap exceeded')
        start=self.pos;end=start+n-1
        req=urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}','User-Agent':'Research-data-verification/1.0'})
        with urllib.request.urlopen(req,timeout=60) as resp:
            if resp.status!=206 or not resp.headers.get('Content-Range','').startswith(f'bytes {start}-'):raise ValueError('Server ignored or changed Range; full download aborted')
            b=resp.read(n+1)
        assert len(b)==n;self.pos+=n;self.bytes+=n;self.requests.append(dict(start=start,end=end,bytes=n));return b

sumatra=json.loads((R/'sumatra_metadata.txt').read_text(encoding='utf-8'))
sumfile=next(f for f in sumatra['files'] if f['key'].startswith('waveforms'))
targets=[('cape','https://gdr.openei.org/files/1823/spectra.zip'),('sumatra',sumfile['links']['self'])]
manifest=[]
for label,url in targets:
    if label=='cape' and (R/'SAMPLE_DOWNLOAD_MANIFEST.json').exists():
        previous=json.loads((R/'SAMPLE_DOWNLOAD_MANIFEST.json').read_text(encoding='utf-8'))
        good=next((a for a in previous if a['dataset']=='cape' and a.get('files')),None)
        if good:manifest.append(good);continue
    try:
        remote=Remote(url)
        with zipfile.ZipFile(remote) as z:
            entries=[i for i in z.infolist() if not i.is_dir()]
            (R/(label+'_archive_index.json')).write_text(json.dumps([dict(name=i.filename,bytes=i.file_size,compressed_bytes=i.compress_size,crc32=f'{i.CRC:08x}') for i in entries],indent=2),encoding='utf-8')
            print(label,'total',remote.size,'members',len(entries),'first',[i.filename for i in entries[:12]],flush=True)
            if label=='cape':selected=[i for i in entries if i.filename.endswith('.npz')][:2]
            else:
                # Fetch records that actually have a published constant-Q result.
                selected=[]
                results=list(csv.DictReader((R/'sumatra/results.csv').open(encoding='utf-8-sig')))
                for row in results:
                    if float(row['alpha'])!=0 or row['component']!='Z' or float(row['snr'])<5:continue
                    prefix=f"waveforms/{row['event']}/{row['event']}.{row['station']}."
                    found=next((i for i in entries if i.filename.startswith(prefix) and '.BHZ.' in i.filename and i.filename.endswith('.sac')),None)
                    if found and found not in selected:selected.append(found)
                    if len(selected)==4:break
                assert len(selected)==4
            files=[]
            for i in selected:
                if i.file_size>8_000_000:continue
                b=z.read(i);dest=R/label/'samples'/Path(i.filename).name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b)
                files.append(dict(member=i.filename,file=str(dest.relative_to(R)),bytes=len(b),crc32=f'{i.CRC:08x}',sha256=hashlib.sha256(b).hexdigest()))
        manifest.append(dict(dataset=label,url=url,archive_bytes=remote.size,transferred_bytes=remote.bytes,requests=remote.requests,files=files))
    except Exception as e:manifest.append(dict(dataset=label,url=url,error=str(e)))
    print(json.dumps(manifest[-1],ensure_ascii=True),flush=True)
(R/'SAMPLE_DOWNLOAD_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
