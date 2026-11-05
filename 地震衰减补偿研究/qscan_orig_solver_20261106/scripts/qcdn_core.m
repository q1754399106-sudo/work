function [r, comp, A] = qcdn_core(dn, dt, freqdom, liter, niter, Q0)
% Operator and solver lines of orig_wang2021/Qcdn.m (lines 43-79) without Madagascar I/O and curvelet thresholding.
% dn: n1 x ntr data; freqdom: Ricker dominant frequency of the inversion wavelet; Q0: constant Q of the operator.
[n1, ntr] = size(dn);
tw=-(n1-1)/2*dt:dt:(n1-1)/2*dt;
wavelet=(1-2*(pi*freqdom.*tw).^2).*exp(-(pi*freqdom.*tw).^2);
wavelet=wavelet./max(abs(wavelet));
MatrixWavelet=zeros(n1,n1);
nn=ceil(n1/2);
for i=1:nn
    MatrixWavelet(1:nn+i-1,i)=wavelet(nn-i+1:end);
end
for i=1:nn-1
    MatrixWavelet(i+1:end,i+nn)=wavelet(1:(end-i));
end
Q = Q0*ones(1,n1);
[~,MatrixAbosrp]=time_absorp_wavelet_matrix(n1,dt,freqdom,Q);
A = MatrixWavelet*MatrixAbosrp;
Param.W=A;
r = zeros(n1,ntr); comp = zeros(n1,ntr);
for tr=1:ntr
    [x,~] = yc_pcg(@mtx_oper,Param,dn(:,tr),zeros(n1,1),liter,niter,false);
    comp(:,tr)=MatrixWavelet*x;
    r(:,tr)=x;
end
comp(isnan(comp))=0;
end
