function run_snr_probe(model_number,snr_db,first,last)
% Addendum A (PLAN_ZH.md): same geology and noise process at a higher Frobenius SNR.
% Estimators: FullCalRefit at c=1, FullCalRefit at the fixed factor 1.2803, and QScan.
P=bv_setup();maxNumCompThreads(P.threads);
s=find(strcmp({P.scenes.model},sprintf('F%d',model_number))&strcmp({P.scenes.condition},'Moderate'));
sc=load_scene(P,s);live=any(sc.D~=0,1);level=find(snr_db==[20 26]);assert(~isempty(level));
out=fullfile(P.root,'results',sprintf('snr_%s_%02ddB_%03d_%03d.mat',sc.model,snr_db,first,last));
assert(~isfile(out),'Partial result exists; refusing to overwrite.');
[n,nx]=size(sc.D0);R=last-first+1;slots=3;fixed=P.grid(8);
M1=zeros(n,nx,slots);M2=M1;REL=zeros(R,slots);replicate=(first:last)';
seed=1e9+1e6*model_number+1e5*(2+level)+replicate;best=zeros(R,1);sigma_mean=best;flags=best;started=tic;
for i=1:R
 D=sc.D0+noise_replica(sc.D0,snr_db,seed(i));sigma=estimate_sigmas(D,sc.cm.t);
 sel=qscan_select(D,sigma,sc.cm,sc.Q0,sc.offset);X=zeros(n,nx,slots);
 [X(:,:,1),a]=full_refit(D,sigma,sc.cm,sc.Q0,sc.offset);[X(:,:,2),b]=full_refit(D,sigma,sc.cm,sc.Q0*fixed,sc.offset);
 if sel.best==8,X(:,:,3)=X(:,:,2);c=b;else,[X(:,:,3),c]=full_refit(D,sigma,sc.cm,sc.Q0*P.grid(sel.best),sc.offset);end
 for f=1:slots,REL(i,f)=norm(X(:,live,f)-sc.T(:,live),'fro')/norm(sc.T(:,live),'fro');end
 M1=M1+X;M2=M2+X.^2;best(i)=sel.best;sigma_mean(i)=mean(sigma);
 flags(i)=a.optimization_flags+b.optimization_flags+c.optimization_flags;
 fprintf('%s %ddB r=%d factor=%.4f rel base %.4f fixed %.4f scan %.4f elapsed %.0fs\n',sc.model,snr_db,replicate(i), ...
  P.grid(best(i)),REL(i,1),REL(i,2),REL(i,3),toc(started));
end
model=sc.model;labels={'c=1.0000','c=1.2803','QScan'};grid=P.grid; %#ok<NASGU>
save(out,'model','snr_db','labels','grid','replicate','seed','best','sigma_mean','flags','REL','M1','M2','-v7');
fprintf('DONE snr probe %s %d dB %.0fs\n',sc.model,snr_db,toc(started));
end
