"""Exploratory (BETA_ZH.md step 1): apparent window t* (field 26 SR lines, semi-synthetic control) vs the curtain model.
Output: analysis/tstar_shape.json."""
import sys, json; sys.path.insert(0,'.')
from env15 import np, lib, semi, gen
sys.path.insert(0,'../../field_baselines_volve_20261029/scripts'); import lib13 as L13
import sq
cm=semi.curtain_cm(); Q=np.asarray(cm['Q'],float); WINS=sq.windows()
def dT(D, c_pred=None):
    sp=[sq.spectra(D,a,b,True) for a,b in WINS]; f=sp[0][0]
    pair=np.logical_and.reduce([s[4] for s in sp]); live=np.logical_and.reduce([s[3] for s in sp]); ip=np.where(pair)[0]; il=np.where(live)[0]
    C=np.array([lib.smooth(s[2][:,ip].mean(1),5) for s in sp]); P=np.array([lib.smooth(s[1][:,il].mean(1),5) for s in sp])
    snr=np.maximum(C,0)/np.maximum(P-np.maximum(C,0),1e-30*P.max(1,keepdims=True))
    fpk=float(f[np.argmax(np.sqrt(np.maximum(C[0],0))*(f>=5))]); hi=min(lib.snr_band(f,snr[-1],fpk),sq.F_CAP)
    sel=(f>=sq.F_LO)&(f<=hi)&np.all(C>0,0); L=0.5*np.log(np.where(C>0,C,np.nan))
    meas=[-np.polyfit(f[sel],(L[k]-L[0])[sel],1)[0]/np.pi for k in range(1,len(WINS))]
    tab=sq.mw_table(cm,Q,f,WINS,grid=np.array([1.0]))[0]   # predicted 0.5 ln power at c=1
    pred1=[-np.polyfit(f[sel],(tab[k]-tab[0])[sel],1)[0]/np.pi for k in range(1,len(WINS))]
    return np.array(meas), np.array(pred1), float(hi)
out={}
SRL=json.load(open(L13.IN/'specslope_q_volve_20261028'/'results'/'real'/'lines.json'))
rows=[]
for ln in sorted(int(k) for k in SRL):
    D=L13.read_inline(ln); m,p1,hi=dT(D); c=SRL[str(ln)]['c']; rows.append((ln,c,m,p1/c,hi))
M=np.array([r[2] for r in rows]); P=np.array([r[3] for r in rows])
print('windows',WINS[1:])
print('FIELD (26 SR lines): measured dT* (ms) mean',np.round(1e3*M.mean(0),2).tolist(),' model at c_SR',np.round(1e3*P.mean(0),2).tolist(),' ratio',np.round((M/P).mean(0),2).tolist())
base=gen.prepare(10227,[10149,10273,10355]); S=[]
for c in (0.7,1.0,1.4):
    for s in (9940,9939):
        g=gen.signal(base,c,s); D=g['x']+g['level']*base['noise'][0]; m,p1,hi=dT(D); S.append((c,m,p1/c))
Ms=np.array([r[1] for r in S]); Ps=np.array([r[2] for r in S])
print('SEMI (6, true c): measured',np.round(1e3*Ms.mean(0),2).tolist(),' model at true c',np.round(1e3*Ps.mean(0),2).tolist(),' ratio',np.round((Ms/Ps).mean(0),2).tolist())
json.dump(dict(windows=WINS,field=[dict(inline=r[0],c_sr=r[1],meas=r[2].tolist(),model_csr=r[3].tolist(),f_hi=r[4]) for r in rows],
               semi=[dict(c=r[0],meas=r[1].tolist(),model=r[2].tolist()) for r in S]),open('../analysis/tstar_shape.json','w'),indent=1)
