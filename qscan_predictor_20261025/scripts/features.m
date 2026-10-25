function features(k)
% Truth-free features Z1-Z6 (PLAN_ZH.md) of model k, all three scenes, replicates 1..R, at the
% QScan choice of each replicate. Writes analysis/features_<model>.csv.
P=qp_setup();mo=P.models(k);[cm,D0]=load_model(P,k);offset=cm.solver_offset;
out=fullfile(P.root,'analysis',sprintf('features_%s.csv',mo.name));assert(~isfile(out),'exists');
fid=fopen(out,'w');fprintf(fid,'model,scene,replicate,seed,best,qscan_factor,Z1,Z2,Z3,Z4,Z5,Z6,seconds\n');
t1=2*cm.h(1)/cm.vp(1);started=tic;
for s=1:numel(P.scenes)
 sc=P.scenes(s);Qused=reshape(cm.Q,1,[])*sc.q_scale;
 for r=1:P.R
  t0=tic;
  if strcmp(mo.family,'G'),seed=P.seed_G(mo.index,sc.level,r);else,seed=P.seed_F(mo.index,sc.level,r);end
  D=D0+noise_replica(D0,sc.snr_db,seed);sig=estimate_sigmas(D,cm.t);sel=qscan_select(D,sig,cm,Qused,offset);
  b=sel.best;J=sel.objective;n=numel(J);
  lo=b-1;hi=b+1;if lo<1,lo=hi;end;if hi>n,hi=lo;end
  Z1=sel.lasso_active(b);Z2=(J(lo)+J(hi)-2*J(b))/J(b);Z3=sel.misfit(b)/J(b);
  S=full_refit(D,sig,cm,Qused*P.grid(b),offset);
  Z4=sqrt(max(norm(D,'fro')^2-cm.nt*sum(sig.^2),0))/norm(S,'fro');
  if b==1||b==n,Z5=0;else,Z5=(J(b-1)-J(b+1))/(2*(J(b-1)+J(b+1)-2*J(b)));end
  w=find(cm.t>=t1);e=round(linspace(0,numel(w),4));E=zeros(1,3);
  for q=1:3,X=D(w(e(q)+1:e(q+1)),:);E(q)=mean(X(:).^2)-mean(sig.^2);end
  Z6=.5*log(max(E(3),1e-12)/E(1));
  fprintf(fid,'%s,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.3f\n',mo.name,sc.tag,r,seed,b,P.grid(b),Z1,Z2,Z3,Z4,Z5,Z6,toc(t0));
  fprintf('%s %s r%d: c=%.4f Z1 %.4f Z2 %.4g Z3 %.4f Z4 %.4f Z5 %+.3f Z6 %+.3f (%.0fs, total %.0fs)\n',mo.name,sc.tag,r,P.grid(b),Z1,Z2,Z3,Z4,Z5,Z6,toc(t0),toc(started));
 end
end
fclose(fid);fprintf('DONE %s\n',mo.name);
end
