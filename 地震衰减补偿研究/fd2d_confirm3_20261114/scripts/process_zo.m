function tr=process_zo(raw,hom,pdir,c0,config,P)
% Zero-offset processing of one 2-D finite-difference trace into the solver's form (plane-wave
% reflection response convolved with the zero-phase Ricker of kennett_trace, 4 ms sampling).
% Applied identically to the attenuated data and to the elastic target.
%  1. direct wave: subtract the source-node trace of the laterally homogeneous cap model;
%  2. signature: divide by the elastic homogeneous trace recorded r0 straight below the source
%     and multiply by the Ricker spectrum and the r0 travel time (removes the source pulse and
%     the 2-D line-source filter in one step);
%  3. 2-D spreading: multiply by sqrt(L(t)/r0), L(t) = (1/c0) * integral of v^2 over two-way
%     time, from the public macro model (config.vp, config.h);
%  4. keep the band below the 4 ms Nyquist and resample to config.t.
% raw, hom, pdir: fd.nt samples at times n*fd.dt (n = 1..fd.nt).
fd=P.fd;dt=fd.dt;dt4=config.dt;r=round(dt4/dt);assert(abs(r*dt-dt4)<1e-12);
N4=2^nextpow2(2*fd.nt/r);N=N4*r;f=(0:N/2)'/(N*dt);
D=fft(raw(:)-hom(:),N);Pd=fft(pdir(:),N);D=D(1:N/2+1);Pd=Pd(1:N/2+1);
src=(2/sqrt(pi))*(f.^2/config.fpeak^3).*exp(-(f/config.fpeak).^2);   % kennett_trace Ricker spectrum
F=src.*exp(-2i*pi*f*fd.r0/c0).*conj(Pd)./(abs(Pd).^2+1e-12*max(abs(Pd).^2));
X=D.*F;                                                           % continuous-FT units of kennett_trace
X4=X(1:N4/2+1);X4(end)=real(X4(end));X4=[X4;conj(X4(end-1:-1:2))];
x=real(ifft(X4))/dt4;x=x(1:config.nt);
% 2-D spreading from the public macro model (velocity as a function of two-way time).
t=config.t;tb=2*cumsum(config.h./config.vp(1:3));
vt=config.vp(1+sum(t>tb(:)',2))';vt=vt(:);
L=cumtrapz(t,vt.^2)/c0;tr=x.*sqrt(L/fd.r0);
end
