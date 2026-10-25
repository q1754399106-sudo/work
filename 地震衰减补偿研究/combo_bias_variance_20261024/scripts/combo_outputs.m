function [G,info]=combo_outputs(P,D,sig,cm,Qused,offset)
% The four estimators of the 2x2 design (plus SoftGroupRefit and QScan+Soft) on one gather,
% statement for statement as round3_20261024/run_combo.m: QScan selection, then for the
% nominal and the scanned Q refit_state -> coherence_scores -> apply_gain. Truth-free.
t0=tic;sig=reshape(sig,1,[]);sel=qscan_select(D,sig,cm,Qused,offset);c=P.grid(sel.best);
G=struct();info=struct('qscan_factor',c,'qscan_best',sel.best,'objective',sel.objective,'optimization_flags',0,'max_normal_residual',0);
for route={'nominal','qscan'}
 if strcmp(route{1},'nominal'),Q=Qused;else,Q=Qused*c;end
 [st,f,nr]=refit_state(D,sig,cm,Q,offset);info.optimization_flags=info.optimization_flags+f;info.max_normal_residual=max(info.max_normal_residual,nr);
 cs=coherence_scores(st.R0,st.supports,st.d,st.live,P.coherence.half_width,P.coherence.dip_periods);
 s=struct('R0',st.R0,'supports',{st.supports},'E',{st.E},'V',{st.V},'groups',{cs.groups},'fallback',{cs.fallback},'W',st.d.W);
 one=struct('lambda_min',1,'lambda_max',1);assert(norm(apply_gain(s,cs.c,one)-st.soft,'fro')<=1e-10*max(norm(st.soft,'fro'),realmin));
 std_=apply_gain(s,cs.c,P.rules.Standard);cc=cell2mat(reshape(cs.c,[],1));
 if strcmp(route{1},'nominal')
  G.FullCalRefit=st.S0;G.Standard=std_;G.SoftGroupRefit=st.soft;info.mean_coherence_nominal=mean(cc);
 else
  G.QScanRefit=st.S0;G.QScan_Standard=std_;G.QScan_Soft=st.soft;info.mean_coherence_qscan=mean(cc);
 end
end
info.seconds=toc(t0);
end

function [st,flags,normal]=refit_state(D,sig,cm,Q,offset)
% refit_state of run_combo.m without its file cache.
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
