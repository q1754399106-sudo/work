function run_trial()
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'source'));
data=jsondecode(fileread(fullfile(root,'inputs','trial.json')));c=struct('dt',data.dt,'nt',data.nt,'nfft',data.nfft,'fref',data.fref);
outdir=fullfile(root,'matlab');assert(~isfolder(outdir));mkdir(outdir);rng(20261023,'twister');started=tic;
allresults=struct([]);total=0;flags=0;maxnormal=0;
for j=1:numel(data.items)
 a=data.items(j);D=a.D(:);N=a.noise(:);sigma=a.sigma;results={};
 for factor=[.75 1 1.25]
  star=a.tstar*factor;G=attenuation_matrix(repmat(star,c.nt,1),c);
  [op,dict]=prepare(G,data.cluster_frequency_hz,c);
  [fixed,refit,pruned,diag]=apply_all(D,op,dict,sigma);
  [fn,rn,pn,dn]=apply_all(N,op,dict,sigma);
  results{end+1}=pack('full_path',factor,star,0,D,[],N,fixed,refit,pruned,fn,rn,pn,G,diag,dn);
  flags=flags+diag.flags+dn.flags;maxnormal=max([maxnormal,diag.normal,dn.normal]);total=total+1;
 end
 for star=[.03 .06 .12]
  G=attenuation_matrix(repmat(star,c.nt,1),c);[op,dict]=prepare(G,data.cluster_frequency_hz,c);
  for nl=[0 1]
   z=randn(c.nt,1);nfft=1024;shape=abs(fft(N,nfft));colored=real(ifft(fft(z,nfft).*shape));colored=colored(1:c.nt);colored=colored-mean(colored);colored=colored/max(std(colored),realmin)*sigma;
   observed=G*D+nl*colored;
   % Pre-existing field noise propagated by G is also present. Keep its sigma
   % in the regularizer; do not exploit the synthetic reference to tune lambda.
   sigma_obs=max(sigma,sqrt((norm(G*N)/sqrt(c.nt))^2+(nl*sigma)^2));
   [fixed,refit,pruned,diag]=apply_all(observed,op,dict,sigma_obs);
   empty=zeros(c.nt,1);
   results{end+1}=pack('controlled_extra',nl,star,nl,observed,D,N,fixed,refit,pruned,empty,empty,empty,G,diag,struct('flags',0,'normal',0));
   flags=flags+diag.flags;maxnormal=max(maxnormal,diag.normal);total=total+1;
  end
 end
 uq_writejson(fullfile(outdir,[a.name '.json']),struct('metadata',rmfield(a,{'D','noise'}),'t',((0:c.nt-1)'*c.dt-4),'observed',D,'background',N,'results',{results}));
 allresults(j).name=a.name;allresults(j).count=numel(results);
 uq_writejson(fullfile(outdir,'progress.json'),struct('records',j,'expected_records',numel(data.items),'cases',total,'elapsed_seconds',toc(started)));
 fprintf('%s %d/%d cases=%d flags=%d elapsed=%.1fs\n',a.station,j,numel(data.items),total,flags,toc(started));
end
uq_writejson(fullfile(outdir,'checks.json'),struct('complete',true,'records',numel(data.items),'cases',total,'optimization_flags',flags,'maximum_normal_residual',maxnormal,'elapsed_seconds',toc(started),'config',data.plan));
end
function [op,d]=prepare(G,freq,c)
I=eye(c.nt);lambda=1e-3;
op=struct('B0',ridge_operator(G,lambda),'G0',G,'BX',{{ridge_operator(I,lambda)}},'GX',{{I}},'lambda_x',lambda,'lambda_fixed',lambda);
scale=sqrt(sum(G.^2,1));B=G./scale;Gram=B'*B;Gram=(Gram+Gram')/2;
d=struct('W',I,'G',G,'B',B,'Gram',Gram,'scale',scale(:),'L',norm(Gram,1),'key',sprintf('%.16g,',[c.dt G(1,1)]),'frequency',freq,'dt',c.dt);
end
function [fixed,refit,pruned,diag]=apply_all(D,op,d,sigma)
[fixed,normal,~]=fixed_v1_apply(D,op);
initial=reflectivity_refit_prepared(D,{d},sigma);out=reflectivity_prune_prepared(D,{d},sigma,initial);
refit=initial.gather;pruned=out.gather;
diag=struct('flags',out.optimization_flags,'normal',max(normal,out.normal_residual),'support_before',out.support_before_pruning,'support_after',out.support_after,'removed',out.removed_coefficients,'source_dictionary','Identity: waveform sparsity field adapter; not Ricker reflectivity');
end
function row=pack(kind,factor,star,nl,D,target,N,fixed,refit,pruned,fn,rn,pn,G,diag,dn)
row=struct('kind',kind,'factor',factor,'tstar',star,'noise_factor',nl,'input',D,'target',target,'background_input',N, ...
 'fixed',fixed,'refit',refit,'optimized',pruned,'fixed_background',fn,'refit_background',rn,'optimized_background',pn, ...
 'fixed_replay',G*fixed,'refit_replay',G*refit,'optimized_replay',G*pruned,'diagnostics',diag,'background_diagnostics',dn);
end
