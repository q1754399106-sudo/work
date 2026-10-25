function run_sweep(s,first,last)
% Stage 2 worker (PLAN_STAGE2_ZH.md): replicates first..last of scene s; r=0 is the noise-free
% gather with the archived sigma. For c=1 and the 13 QScan grid factors, the full pipeline
% (refit_state -> coherence_scores -> apply_gain Standard) at Q = Qused*c, plus FullCalRefit at the
% same Q. r>=1 is checked bit for bit (as single) against the stage-1 cache. Truth is not read.
P=cb_setup();maxNumCompThreads(P.threads);sc=load_scene(P,s);assert(first>=0&&last<=20);
dir_=fullfile(P.root,'cache','sweep',sc.name);if ~isfolder(dir_),mkdir(dir_);end
factors=[1 P.grid];nf=numel(factors);[n,nx]=size(sc.D0);started=tic;
for r=first:last
 out=fullfile(dir_,sprintf('r%03d.mat',r));if isfile(out),continue;end
 t0=tic;
 if r==0
  in=load(fullfile(P.bae,'inputs','independent',[sc.model '.mat']),'obs');cond='Moderate';if sc.snr_db==8,cond='Strong';end
  D=sc.D0;sig=reshape(in.obs(strcmp({in.obs.name},cond)).sigma,1,[]);seed=NaN;
 else
  seed=sc.seed_base+r;D=sc.D0+noise_replica(sc.D0,sc.snr_db,seed);sig=estimate_sigmas(D,sc.cm.t);
 end
 sel=qscan_select(D,sig,sc.cm,sc.Qused,sc.offset);STD=zeros(n,nx,nf,'single');FCR=STD;flags=0;normal=0;
 for f=1:nf
  [st,fl,nr]=refit_state(D,sig,sc.cm,sc.Qused*factors(f),sc.offset);flags=flags+fl;normal=max(normal,nr);
  cs=coherence_scores(st.R0,st.supports,st.d,st.live,P.coherence.half_width,P.coherence.dip_periods);
  g=struct('R0',st.R0,'supports',{st.supports},'E',{st.E},'V',{st.V},'groups',{cs.groups},'fallback',{cs.fallback},'W',st.d.W);
  STD(:,:,f)=single(apply_gain(g,cs.c,P.rules.Standard));FCR(:,:,f)=single(st.S0);
 end
 check=NaN;
 if r>=1
  L=load(fullfile(P.root,'cache','outputs',sc.name,sprintf('r%03d.mat',r)),'G','info');assert(L.info.qscan_best==sel.best);
  check=isequal(STD(:,:,1),L.G.Standard)&&isequal(FCR(:,:,1),L.G.FullCalRefit)&&isequal(STD(:,:,1+sel.best),L.G.QScan_Standard);
  assert(check,'Sweep does not reproduce the stage-1 outputs.');
 end
 info=struct('factors',factors,'objective',sel.objective,'best',sel.best,'optimization_flags',flags,'max_normal_residual',normal, ...
  'stage1_identical',check,'seconds',toc(t0));
 save([out '.part'],'STD','FCR','info','seed','sig','-v7');movefile([out '.part'],out);
 fprintf('%s r=%d best %.4f flags %d stage1 %d (%.0fs, total %.0fs)\n',sc.name,r,P.grid(sel.best),flags,check,info.seconds,toc(started));
end
fprintf('DONE sweep %s %d-%d %.0fs\n',sc.name,first,last,toc(started));
end

function [st,flags,normal]=refit_state(D,sig,cm,Q,offset)
% refit_state of round3_20261024/run_combo.m without its file cache (as in combo_outputs.m).
c=cm;c.Q=Q;d=prune_dictionary(c,offset);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;dicts=repmat({d},1,11);nx=size(D,2);
st=struct('R0',zeros(size(D)),'S0',zeros(size(D)),'soft',zeros(size(D)),'supports',{cell(1,nx)},'E',{cell(1,nx)},'V',{cell(1,nx)},'live',any(D~=0,1));
flags=0;normal=0;
for b=1:ceil(nx/11)
 ids=(b-1)*11+1:min(b*11,nx);ini=reflectivity_refit_prepared(D(:,ids),dicts(1:numel(ids)),sig(ids));so=reflectivity_soft_prepared(D(:,ids),dicts(1:numel(ids)),sig(ids),ini);
 st.R0(:,ids)=ini.reflectivity;st.S0(:,ids)=ini.gather;st.soft(:,ids)=so.gather;st.supports(ids)=ini.support;st.E(ids)=so.group_energy;st.V(ids)=so.group_noise_energy;
 flags=flags+ini.optimization_flags;normal=max([normal,ini.normal_residual,so.initial_normal_residual]);
end
st.d=d;
end
