function G = attenuation_matrix(a,c)
% Column j is the response of an impulse at time t_j.
% MATLAB FFT convention: forward exp(-i*2*pi*f*t).
% Positive-frequency H=exp(-pi*f*a + i*2*f*log(f/fref)*a).
% Negative frequencies use conjugate symmetry; DC=1, Nyquist real.
N=c.nfft; f=(0:N/2)'/(N*c.dt);
b=-pi*f+1i*2*f.*log(max(f,realmin)/c.fref); b(1)=0;
H=exp(b*a(:)'); H(end,:)=real(H(end,:));
H=[H; conj(H(end-1:-1:2,:))];
E=fft(eye(c.nt),N,1);
Z=ifft(H.*E,[],1); G=real(Z(1:c.nt,:));
end
