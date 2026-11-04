import numpy as np
N=256; dt=0.002; Q=40
t=np.arange(N)*dt
f=np.fft.rfftfreq(N,dt)
def atten_matrix(Q):
    A=np.zeros((N,N))
    for j,tau in enumerate(t):
        H=np.exp(-np.pi*f*tau/Q)*np.exp(-1j*2*np.pi*f*tau)
        A[:,j]=np.fft.irfft(H,n=N)
    return A
def shift(g):
    return np.roll(np.eye(N),g,axis=0)
def erank(M,thr):
    s=np.linalg.svd(M,compute_uv=False); return int((s>thr*s.max()).sum()), s
A=atten_matrix(Q)
# stationary circulant conv with Ricker
fp=40; tw=(np.arange(N)-N//2)*dt
w=(1-2*(np.pi*fp*tw)**2)*np.exp(-(np.pi*fp*tw)**2); w=np.roll(w,-N//2)
C=np.stack([np.roll(w,j) for j in range(N)],axis=1)
shifts=list(range(0,N,8))
for thr in [1e-2,1e-3]:
    rA,_=erank(A,thr)
    S=np.vstack([A@shift(g) for g in shifts])/np.sqrt(len(shifts))
    rS,_=erank(S,thr)
    rC,_=erank(C,thr)
    SC=np.vstack([C@shift(g) for g in shifts])/np.sqrt(len(shifts))
    rSC,_=erank(SC,thr)
    print(f"thr={thr}: atten A rank {rA} -> stacked {rS} / {N};  stationary conv rank {rC} -> stacked {rSC}")
# time-frequency check: recoverable bandwidth at t=0 equals? compare with A restricted
# Non-circular shifts (zero fill) limited to +-80 samples
def shift_nc(g):
    M=np.zeros((N,N))
    for i in range(N):
        j=i-g
        if 0<=j<N: M[i,j]=1
    return M
for thr in [1e-2,1e-3]:
    S=np.vstack([A@shift_nc(g) for g in range(-80,81,8)])
    print("noncircular +-80:",thr, erank(S,thr)[0])
