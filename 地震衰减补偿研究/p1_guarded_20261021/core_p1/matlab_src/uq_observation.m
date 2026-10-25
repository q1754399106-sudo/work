function [D,sigma]=uq_observation(scene,c,trueModels,map,seed,snr,stress)
% SYNTHESIS ONLY. No inverse method may call this function.
if nargin<7, stress=[]; end
signal=scene.truth;
if ~isempty(stress) && strcmp(stress.name,'StressNonRickerPhase')
 z=fft(signal,c.nfft,1); rad=stress.phase_deg*pi/180;
 ph=[1;exp(-1i*rad)*ones(c.nfft/2-1,1);1;exp(1i*rad)*ones(c.nfft/2-1,1)];
 signal=real(ifft(z.*ph,[],1)); signal=signal(1:c.nt,:);
end
clean=zeros(size(signal)); for j=1:11, clean(:,j)=trueModels{map(j)}.G*signal(:,j); end
if ~isempty(stress) && strcmp(stress.name,'StressResidualNMO')
 f=[0:c.nfft/2,-c.nfft/2+1:-1]'/(c.nfft*c.dt);
 for j=1:11
  shift=stress.max_nmo_shift_s*((j-1)/10)^2;
  z=real(ifft(fft(clean(:,j),c.nfft).*exp(-2i*pi*f*shift))); clean(:,j)=z(1:c.nt);
 end
end
coherent=zeros(size(clean));
if strcmp(scene.name,'CoherentNoiseAVO')
 for j=1:11, z=pi*35*(c.t-(mean(c.event_t(1:2))+.004*(j-1))); coherent(:,j)=(1-2*z.^2).*exp(-z.^2); end
 coherent=.25*norm(clean,'fro')*coherent/norm(coherent,'fro');
end
rng(seed,'twister'); noise=randn(size(clean));
if ~isempty(stress) && strcmp(stress.name,'StressColoredNoise')
 rho=stress.ar1_rho;
 for k=2:size(noise,1), noise(k,:)=rho*noise(k-1,:)+sqrt(1-rho^2)*noise(k,:); end
end
noise=noise./sqrt(mean(noise.^2,1));
sigma=sqrt(mean(clean.^2,1))*10^(-snr/20); D=clean+noise.*sigma+coherent;
end
