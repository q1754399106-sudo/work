function S=wang_inverse_q(D,cm,offset,gain_limit_db)
% Stabilized inverse Q filter (Wang 2006, Geophysics): unconditionally stable
% phase operator; amplitude operator (beta+s2)/(beta^2+s2), s2=exp(-(0.23*G+1.63)).
% Uses the same absorption integral and dispersion law as attenuation_matrix.
t=cm.t(:);n=numel(t);N=cm.nfft;[a,b]=ray_absorption(t,offset,cm);tau=a+b;
f=(0:N/2)'/(N*cm.dt);s2=exp(-(.23*gain_limit_db+1.63));
beta=exp(-pi*f*tau');phase=2*f.*log(max(f,realmin)/cm.fref)*tau';phase(1,:)=0;
H=((beta+s2)./(beta.^2+s2)).*exp(-1i*phase);H(end,:)=real(H(end,:));H=[H;conj(H(end-1:-1:2,:))];
m=(0:N-1)';P=(H.*exp(1i*2*pi*m*(0:n-1)/N)).'/N;
S=real(P*fft(D,N,1));
end
