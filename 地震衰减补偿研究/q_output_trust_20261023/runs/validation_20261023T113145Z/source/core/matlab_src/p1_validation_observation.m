function [D,sigma]=p1_validation_observation(scene,c,trueModels,map,seed,su)
if any(strcmp(su.kind,{'normal','velocity'}))
 [D,sigma]=uq_observation(scene,c,trueModels,map,seed,su.snr_db); return;
end
[n,nx]=size(scene.truth); clean=zeros(n,nx); N=c.nfft;
for j=1:nx
 if strcmp(su.kind,'phase')
  rad=su.phase_deg*pi/180;
  phase=[1;exp(-1i*rad)*ones(N/2-1,1);1;exp(1i*rad)*ones(N/2-1,1)];
  wr=real(ifft(fft(scene.truth(:,j),N).*phase)); clean(:,j)=trueModels{map(j)}.G*wr(1:n);
 else
  raw=trueModels{map(j)}.G*scene.truth(:,j);
  if strcmp(su.kind,'nmo')
   f=(0:N/2)'/(N*c.dt); wf=fft(raw,N); shift=su.max_nmo_shift_s*((j-1)/(nx-1))^2;
   half=wf(1:N/2+1).*exp(-2i*pi*f*shift);
   padded=real(ifft([half;conj(half(end-1:-1:2))])); raw=padded(1:n);
  end
  clean(:,j)=raw;
 end
end
coherent=zeros(n,nx);
if strcmp(scene.name,'CoherentNoiseAVO')
 for j=1:nx, z=pi*35*(c.t-(mean(c.event_t(1:2))+.004*(j-1))); coherent(:,j)=(1-2*z.^2).*exp(-z.^2); end
 coherent=.25*norm(clean,'fro')*coherent/norm(coherent,'fro');
end
rng(seed,'twister'); noise=randn(n,nx);
if strcmp(su.kind,'colored')
 rho=su.ar1_rho;
 for t=2:n, noise(t,:)=rho*noise(t-1,:)+sqrt(1-rho^2)*noise(t,:); end
end
noise=noise./sqrt(mean(noise.^2,1)); sigma=sqrt(mean(clean.^2,1))*10^(-su.snr_db/20);
D=clean+noise.*sigma+coherent;
end
