"""Task-local I/O; immutable Git/LFS inputs, single BLAS thread."""
import os
for _k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_k]='1'
from pathlib import Path
import sys, subprocess, hashlib, json
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parent
sys.path.insert(0,str(REPO/'novelty_evidence_20261025'/'cache'/'deps'))
import numpy as np
from scipy.io import loadmat
from scipy import sparse

def gitbytes(path):
    raw=subprocess.run(['git','-C',str(REPO),'show',f'198f04b76d54c5ca77cd7b59df3a5653d6001a3e:{path}'],capture_output=True,check=True).stdout
    if raw.startswith(b'version https://git-lfs.github.com/spec/v1'):
        f=dict(x.split(' ',1) for x in raw.decode().strip().splitlines())
        oid=f['oid'].split(':')[1]
        raw=(REPO/'.git/lfs/objects'/oid[:2]/oid[2:4]/oid).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==oid and len(raw)==int(f['size'])
    return raw

def archived_mat(path):
    p=ROOT/'cache'/'inputs'/path
    if not p.exists():
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(gitbytes(path))
    return loadmat(p,simplify_cells=True)

def dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)

def vec(x,dtype=float):
    return np.asarray(x,dtype=dtype).reshape(-1)

def writejson(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

def dictionary(cm):
    t=vec(cm['t']);n=len(t);dt=float(cm['dt']);off=float(cm.get('solver_offset',450));q=vec(cm['Q']);v=vec(cm['vp']);h=vec(cm['h'])
    bounds=np.r_[0,np.cumsum(2*h/v[:3]),np.inf]
    tau=np.maximum(0,np.minimum(t[:,None],bounds[1:])-bounds[:-1])
    a=np.sum(tau/q,axis=1)
    if off:
        for k in range(n):
            active=tau[k]>0
            if not active.any(): a[k]=off/(v[0]*q[0]);continue
            vv=v[active];tt=tau[k,active];qq=q[active];lo=0.;hi=(1-1e-12)/vv.max()
            for _ in range(70):
                p=(lo+hi)/2
                if np.sum(tt*vv**2*p/np.sqrt(1-(p*vv)**2))<off:lo=p
                else:hi=p
            p=(lo+hi)/2;a[k]=np.sum(tt/qq/np.sqrt(1-(p*vv)**2))
    N=int(cm['nfft']);f=np.arange(N//2+1)/(N*dt)
    b=-np.pi*f+2j*f*np.log(np.maximum(f,np.finfo(float).tiny)/float(cm['fref']));b[0]=0
    G=np.empty((n,n))
    # Column streaming avoids nfft-by-nt complex temporaries (>100MB for V5).
    for k in range(n):
        H=np.exp(b*a[k]);H[-1]=H[-1].real
        phase=np.exp(-2j*np.pi*np.arange(N//2+1)*k/N)
        G[:,k]=np.fft.irfft(H*phase,n=N)[:n]
    z=np.pi*float(cm['fpeak'])*(t[:,None]-t[None,:])
    W=(1-2*z*z)*np.exp(-z*z)
    B=G@W;scale=np.linalg.norm(B,axis=0);B/=scale
    return B,W,W/scale,scale

def load_scene(name,noise):
    a=archived_mat(f'round3_20261024/dev/{name}.mat')
    tr=archived_mat(f'round3_20261024/dev/{name}_truth.mat')
    cm=a['config'];obs=a['obs'];obs=obs if isinstance(obs,list) else [obs]
    o=next(x for x in obs if x['name']==noise)
    raw=ROOT/f'inputs/prepared/{name}_Q+0_{noise}_nominal.mat'
    if not raw.exists():
        raw.parent.mkdir(parents=True,exist_ok=True)
        raw.write_bytes(prepared_bytes(name,noise))
    if not raw.is_file():raise FileNotFoundError(raw)
    st=loadmat(raw,simplify_cells=True)['st']
    return cm,np.asarray(o['D'],float),vec(o['sigma']),st,tr


def prepared_bytes(name,noise):
    return pinned_bytes(f"novelty_evidence_20261025/inputs/prepared/{name}_Q+0_{noise}_nominal.mat")

def pinned_bytes(path):
    from bootstrap import gitbytes as read_pinned
    return read_pinned(path)
