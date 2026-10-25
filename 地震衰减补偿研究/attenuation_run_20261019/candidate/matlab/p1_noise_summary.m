function out=p1_noise_summary(D,t,o,seed)
% Joint circular block bootstrap preserves within-window temporal/trace pairing.
% An upper empirical quantile is a conservative plug-in, not a confidence limit.
ix=t>=o.noise_window(1) & t<=o.noise_window(2);
Z=D(ix,:); assert(size(Z,1)>=32,'Noise window needs at least 32 samples');
Z=Z-median(Z,1); n=size(Z,1); nx=size(D,2); c=.6744897501960817;
sigma=median(abs(Z),1)/c;
floorSigma=max(sqrt(mean(D.^2,1))*1e-12,sqrt(realmin));
sigma=max(sigma,floorSigma); rho=zeros(1,nx);
for j=1:nx
 z=Z(:,j); z=z-mean(z); rho(j)=sum(z(1:end-1).*z(2:end))/max(sum(z.^2),realmin);
end
rho=max(-o.rho_limit,min(o.rho_limit,rho));
stream=RandStream('mt19937ar','Seed',seed); samples=zeros(o.noise_bootstrap_replicates,nx);
for k=1:o.noise_bootstrap_replicates
 ids=p1_block_indices(n,n,o.block_length,stream); zz=Z(ids,:);
 samples(k,:)=max(median(abs(zz-median(zz,1)),1)/c,floorSigma).^2;
end
sorted=sort(samples,1); q=1+(size(sorted,1)-1)*o.noise_quantile;
lo=floor(q); hi=ceil(q); v=(hi-q)*sorted(lo,:)+(q-lo)*sorted(hi,:);
if lo==hi, v=sorted(lo,:); end
% Do not call a below-point-estimate bootstrap fluctuation "conservative".
v=max(v,sigma.^2);
out=struct('sigma',sigma,'variance_reference',sigma.^2,'variance_upper',v, ...
 'bootstrap_variances',samples,'rho',rho,'noise_samples',n,'window',o.noise_window, ...
 'seed',seed,'quantile',o.noise_quantile,'block_length',o.block_length, ...
 'upper_quantile_is_confidence_interval',false,'centered_window',Z);
end
