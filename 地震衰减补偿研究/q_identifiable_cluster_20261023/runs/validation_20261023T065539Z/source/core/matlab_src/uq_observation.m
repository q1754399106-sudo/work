function [D,sigma]=uq_observation(scene,c,trueModels,map,seed,snr)
% SYNTHESIS ONLY. No inverse method may call this function.
clean=zeros(size(scene.truth)); for j=1:11, clean(:,j)=trueModels{map(j)}.G*scene.truth(:,j); end
coherent=zeros(size(clean));
if strcmp(scene.name,'CoherentNoiseAVO')
 for j=1:11, z=pi*35*(c.t-(mean(c.event_t(1:2))+.004*(j-1))); coherent(:,j)=(1-2*z.^2).*exp(-z.^2); end
 coherent=.25*norm(clean,'fro')*coherent/norm(coherent,'fro');
end
rng(seed,'twister'); noise=randn(size(clean)); noise=noise./sqrt(mean(noise.^2,1));
sigma=sqrt(mean(clean.^2,1))*10^(-snr/20); D=clean+noise.*sigma+coherent;
end
