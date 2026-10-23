import numpy as np
def soft_fk(D,dt,dx):
 """Declared baseline, no reference truth: low-frequency slow-velocity cone."""
 D=np.asarray(D,float);nt,nx=D.shape;N=2**int(np.ceil(np.log2(2*nt)));M=2**int(np.ceil(np.log2(2*nx)))
 f=np.fft.rfftfreq(N,dt);k=np.fft.fftfreq(M,dx);F=np.fft.rfft(D,n=N,axis=0);power=np.mean(abs(F)**2,axis=1);q=(f>=2)&(f<=min(.45/dt,120));cum=np.cumsum(np.where(q,power,0));f90=float(f[np.searchsorted(cum,.9*cum[-1])]) if cum[-1]>0 else 30.
 pad=max(1,int(round(.08*nx)));taper=np.ones(nx);ramp=.5*(1-np.cos(np.linspace(0,np.pi,pad)));taper[:pad]=ramp;taper[-pad:]=ramp[::-1]
 S=np.fft.fft(np.fft.rfft(D*taper[None,:],n=N,axis=0),n=M,axis=1)
 p=np.divide(abs(k)[None,:],f[:,None],out=np.zeros((len(f),M)),where=f[:,None]>0)
 slope=np.clip((p-1/2500)/(1/1200-1/2500),0,1);slope=.5-.5*np.cos(np.pi*slope)
 freq=np.clip((f-.4*f90)/(.4*f90),0,1);freq=.5+.5*np.cos(np.pi*freq)
 mask=.6*slope*freq[:,None];mask[0,:]=0
 noise=np.fft.irfft(np.fft.ifft(S*mask,axis=1)[:,:nx],n=N,axis=0)[:nt]
 live=np.sum(D*D,axis=0)>0;noise[:,~live]=0;out=D-noise
 return out,noise,dict(method='Soft f-k baseline',full_suppression_below_velocity_m_s=1200,transition_up_to_velocity_m_s=2500,frequency_pass_transition_relative_f90=[.4,.8],maximum_subtraction_fraction=.6,observed_f90_hz=f90,removed_norm_fraction=float(np.linalg.norm(noise)/np.linalg.norm(D)),note='Suppressed spectral cone can contain true signal. Not a noise-free reference and not proof of SNR improvement.')
