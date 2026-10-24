function run_stage0()
% Stage 0 (DESIGN_ZH.md), before freezing: (1) the round-3 archived F1/F2 gathers at Q0 and
% Q-30% reproduce combo_metrics.csv; (2) replicates 1-2 of S1-S4 reproduce FullCalRefit and
% QScan of qscan_bias_variance replicates.csv, and the refit_state path equals full_refit;
% (3) timing and MATLAB memory per replicate for every scene type.
P=cb_setup();maxNumCompThreads(P.threads);A=fullfile(P.root,'analysis');if ~isfolder(A),mkdir(A);end
ref=readtable(fullfile(P.r3,'analysis','combo_metrics.csv'),'TextType','string');rows={};started=tic;
for model={'F1','F2'}
 r3=load(fullfile(P.r3,'dev',[model{1} '.mat']));r3t=load(fullfile(P.r3,'dev',[model{1} '_truth.mat']));
 bae=load(fullfile(P.bae,'inputs','independent',[model{1} '.mat']));
 cm=r3.config;cm.t=cm.t(:);Q0=reshape(cm.Q,1,[]);offset=450;if isfield(cm,'solver_offset'),offset=cm.solver_offset;end
 for k=1:numel(r3.obs)
  o=r3.obs(k);b=bae.obs(strcmp({bae.obs.name},o.name));same=isequal(o.D,b.D)&&isequal(o.sigma,b.sigma);
  T=r3t.tr(strcmp({r3t.tr.name},o.name)).target;
  for e=[0 -30]
   [G,info]=combo_outputs(P,o.D,reshape(o.sigma,1,[]),cm,Q0*(1+e/100),offset);
   for m=P.methods
    x=norm(G.(m{1})-T,'fro')/norm(T,'fro');
    a=ref.overall(ref.set==model{1}&ref.condition==o.name&ref.Q_percent==e&ref.method==m{1});
    f=ref.factor(ref.set==model{1}&ref.condition==o.name&ref.Q_percent==e&ref.method==m{1});
    rows(end+1,:)={'archived_round3',model{1},o.name,e,0,m{1},info.qscan_factor,f,x,a,abs(x-a)/a,same,info.seconds,info.optimization_flags}; %#ok<AGROW>
   end
   fprintf('archived %s %s Q%+d: factor %.4f (%.0fs, total %.0fs)\n',model{1},o.name,e,info.qscan_factor,info.seconds,toc(started));
  end
 end
end
rep=readtable(fullfile(P.bv,'analysis','replicates.csv'),'TextType','string');mem=zeros(0,1);
for s=[1:4 5 9]
 sc=load_scene(P,s);ET=sum(sc.T(:).^2);cond="Moderate";if sc.snr_db==8,cond="Strong";end
 for r=1:2
  if s>4&&r>1,continue;end
  D=sc.D0+noise_replica(sc.D0,sc.snr_db,sc.seed_base+r);sig=estimate_sigmas(D,sc.cm.t);
  [G,info]=combo_outputs(P,D,sig,sc.cm,sc.Qused,sc.offset);
  base=full_refit(D,sig,sc.cm,sc.Qused,sc.offset);scan=full_refit(D,sig,sc.cm,sc.Qused*info.qscan_factor,sc.offset);
  path_ok=isequal(base,G.FullCalRefit)&&isequal(scan,G.QScanRefit);
  [u,~]=memory;mem(end+1,1)=u.MemUsedMATLAB/2^20; %#ok<AGROW>
  for m=P.methods
   x=sqrt(sum((G.(m{1})-sc.T).^2,'all')/ET);a=NaN;f=NaN;
   if s<=4
    q=rep(rep.model==sc.model&rep.condition==cond&rep.replicate==r,:);f=q.qscan_factor;
    if strcmp(m{1},'FullCalRefit'),a=q.rel_c1p0000;elseif strcmp(m{1},'QScanRefit'),a=q.rel_QScan;end
   end
   rows(end+1,:)={'replicate',sc.model,sc.name,round(100*sc.q_scale)-100,r,m{1},info.qscan_factor,f,x,a,abs(x-a)/a,path_ok,info.seconds,info.optimization_flags}; %#ok<AGROW>
  end
  fprintf('%s r=%d: factor %.4f path %d mem %.0f MB (%.0fs, total %.0fs)\n',sc.name,r,info.qscan_factor,path_ok,mem(end),info.seconds,toc(started));
 end
end
writetable(cell2table(rows,'VariableNames',{'check','model','condition','Q_percent','replicate','method','factor','reference_factor', ...
 'rel_error','reference_rel_error','relative_difference','inputs_or_path_identical','seconds','optimization_flags'}),fullfile(A,'stage0_checks.csv'));
fprintf('DONE stage 0 %.0fs, max MATLAB memory %.0f MB\n',toc(started),max(mem));
end
