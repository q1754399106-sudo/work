function run_fd(m,s,first,last)
% Worker (PLAN_ZH.md): replicates first..last of model m, scene s. r=0 is the noise-free gather
% with sigma estimated on replicate 1 of the same scene. For c=1 and the 13 QScan grid factors,
% the full pipeline (refit_state -> coherence_scores -> apply_gain Standard) at Q = Qused*c, plus
% FullCalRefit at the same Q; QScan's objective and choice on the same gather. Truth is not read.
% Statement for statement as qscan_offset_20261025/scripts/run_offset.m, with fd_setup in place
% of os_setup.
P=fd_setup();maxNumCompThreads(P.threads);sp=P.models(m);sc=P.scenes(s);assert(first>=0&&last<=P.R);
M=load(fullfile(P.root,'models',[sp.name '.mat']));cm=M.config;cm.t=cm.t(:);D0=M.clean;
Qused=reshape(cm.Q,1,[])*sc.q_scale;offset=cm.solver_offset;
name=[sp.name '_' sc.tag];dir_=fullfile(P.root,'cache','sweep',name);if ~isfolder(dir_),mkdir(dir_);end
factors=P.factors;nf=numel(factors);[n,nx]=size(D0);started=tic;
for r=first:last
 out=fullfile(dir_,sprintf('r%03d.mat',r));if isfile(out),continue;end
 t0=tic;
 if r==0
  D=D0;seed=NaN;sig=estimate_sigmas(D0+noise_replica(D0,sc.snr_db,P.seed(m,sc.level,1)),cm.t);
 else
  seed=P.seed(m,sc.level,r);D=D0+noise_replica(D0,sc.snr_db,seed);sig=estimate_sigmas(D,cm.t);
 end
 sel=qscan_select(D,sig,cm,Qused,offset);STD=zeros(n,nx,nf,'single');FCR=STD;flags=0;normal=0;
 for f=1:nf
  [st,fl,nr]=refit_state(D,sig,cm,Qused*factors(f),offset);flags=flags+fl;normal=max(normal,nr);
  cs=coherence_scores(st.R0,st.supports,st.d,st.live,P.coherence.half_width,P.coherence.dip_periods);
  g=struct('R0',st.R0,'supports',{st.supports},'E',{st.E},'V',{st.V},'groups',{cs.groups},'fallback',{cs.fallback},'W',st.d.W);
  STD(:,:,f)=single(apply_gain(g,cs.c,P.rules.Standard));FCR(:,:,f)=single(st.S0);
 end
 info=struct('factors',factors,'objective',sel.objective,'best',sel.best,'optimization_flags',flags,'max_normal_residual',normal,'seconds',toc(t0));
 save([out '.part'],'STD','FCR','info','seed','sig','-v7');movefile([out '.part'],out);
 fprintf('%s r=%d best %.4f flags %d (%.0fs, total %.0fs)\n',name,r,P.grid(sel.best),flags,info.seconds,toc(started));
end
fprintf('DONE %s %d-%d %.0fs\n',name,first,last,toc(started));
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
