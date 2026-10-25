function run_replicates(s,first,last)
% Stage 2 worker: replicates first..last of scene s. Each replicate draws a new noise field
% from the archived generator, re-estimates sigma from the noisy gather, and runs QScan and
% FullCalRefit at every factor on that same gather (paired design). Truth is used only to score.
P=bv_setup();maxNumCompThreads(P.threads);sc=load_scene(P,s);live=any(sc.D~=0,1);
out=fullfile(P.root,'results',sprintf('rep_%s_%s_%03d_%03d.mat',sc.model,sc.condition,first,last));
assert(~isfile(out),'Partial result exists; refusing to overwrite.');
[n,nx]=size(sc.D0);nf=numel(P.factors);slots=nf+1;R=last-first+1;   % slot nf+1 = QScan
M1=zeros(n,nx,slots);M2=M1;C=zeros(n,nx);
ids=unique([1:3:99,23,78]);edges=.30:.05:1.60;nw=numel(edges)-1;
REL=zeros(R,slots);SUP=REL;OBJ=zeros(R,numel(P.grid));MIS=OBJ;PEN=OBJ;ACT=OBJ;WM=zeros(R,numel(ids)*nw,2);
replicate=(first:last)';seed=sc.seed_base+replicate;best=zeros(R,1);sigma_mean=best;sigma_cv=best;flags=best;normal=best;
started=tic;
for i=1:R
 D=sc.D0+noise_replica(sc.D0,sc.snr_db,seed(i));sigma=estimate_sigmas(D,sc.cm.t);
 sel=qscan_select(D,sigma,sc.cm,sc.Q0,sc.offset);X=zeros(n,nx,slots);
 for f=1:nf
  [X(:,:,f),info]=full_refit(D,sigma,sc.cm,sc.Q0*P.factors(f),sc.offset);
  SUP(i,f)=mean(info.support_size);flags(i)=flags(i)+info.optimization_flags;normal(i)=max(normal(i),info.max_normal_residual);
 end
 X(:,:,slots)=X(:,:,1+sel.best);SUP(i,slots)=SUP(i,1+sel.best);
 for f=1:slots,REL(i,f)=norm(X(:,live,f)-sc.T(:,live),'fro')/norm(sc.T(:,live),'fro');end
 M1=M1+X;M2=M2+X.^2;C=C+X(:,:,1).*X(:,:,slots);
 w=0;
 for j=ids
  for k=1:nw
   w=w+1;m=sc.cm.t>=edges(k)-1e-9&sc.cm.t<edges(k+1)-1e-9;
   WM(i,w,1)=mean((X(m,j,1)-sc.T(m,j)).^2);WM(i,w,2)=mean((X(m,j,slots)-sc.T(m,j)).^2);
  end
 end
 best(i)=sel.best;OBJ(i,:)=sel.objective;MIS(i,:)=sel.misfit;PEN(i,:)=sel.penalty;ACT(i,:)=sel.lasso_active;
 sigma_mean(i)=mean(sigma);sigma_cv(i)=std(sigma)/mean(sigma);
 fprintf('%s %s r=%d factor=%.4f rel base %.4f scan %.4f flags %d elapsed %.0fs\n',sc.model,sc.condition, ...
  replicate(i),P.grid(best(i)),REL(i,1),REL(i,slots),flags(i),toc(started));
end
model=sc.model;condition=sc.condition;factors=P.factors;grid=P.grid; %#ok<NASGU>
save(out,'model','condition','factors','grid','replicate','seed','best','sigma_mean','sigma_cv','flags','normal', ...
 'REL','SUP','OBJ','MIS','PEN','ACT','WM','M1','M2','C','ids','edges','-v7');
fprintf('DONE stage 2 %s %s %d-%d %.0fs\n',sc.model,sc.condition,first,last,toc(started));
end
