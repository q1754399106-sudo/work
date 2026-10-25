function diag_threads()
% Does the Standard mismatch to round 3 depend on the thread count? F2/Moderate Q-30 and
% F2/Strong Q-30 at 1 and 4 threads; also how much the coherence scores move between them.
P=cb_setup();
r3=load(fullfile(P.r3,'dev','F2.mat'));r3t=load(fullfile(P.r3,'dev','F2_truth.mat'));
ref=readtable(fullfile(P.r3,'analysis','combo_metrics.csv'),'TextType','string');
cm=r3.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);offset=cm.solver_offset;
for k=1:2
 o=r3.obs(k);T=r3t.tr(k).target;
 for th=[1 4]
  maxNumCompThreads(th);[G,info]=combo_outputs(P,o.D,reshape(o.sigma,1,[]),cm,Q0*(1+(-30)/100),offset);
  for m={'Standard','QScan_Standard'}
   x=norm(G.(m{1})-T,'fro')/norm(T,'fro');a=ref.overall(ref.set=="F2"&ref.condition==o.name&ref.Q_percent==-30&ref.method==m{1});
   fprintf('%s threads %d %s: %.9f ref %.9f reldiff %.2e coh %.9f\n',o.name,th,m{1},x,a,abs(x-a)/a,info.mean_coherence_nominal);
  end
  if th==1,S1=G.Standard;else,fprintf('  1 vs 4 threads Standard rel diff %.2e\n',norm(G.Standard-S1,'fro')/norm(S1,'fro'));end
 end
end
end
