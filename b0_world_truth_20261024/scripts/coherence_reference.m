function [XC,x1,info]=coherence_reference(P,D,sig,cm,Q,offset)
% Coherence-shrunk output ("Standard", lambda 0.5..16) on the FullCalRefit groups, exactly as
% round3_20261024/run_combo.m: refit_state -> coherence_scores -> apply_gain. Truth-free.
c=cm;c.Q=Q;d=prune_dictionary(c,offset);Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;
nx=size(D,2);R0=zeros(size(D));x1=R0;soft=R0;supports=cell(1,nx);E=supports;V=supports;sig=reshape(sig,1,[]);
for first=1:11:nx
 ids=first:min(first+10,nx);dicts=repmat({d},1,numel(ids));
 ini=reflectivity_refit_prepared(D(:,ids),dicts,sig(ids));so=reflectivity_soft_prepared(D(:,ids),dicts,sig(ids),ini);
 R0(:,ids)=ini.reflectivity;x1(:,ids)=ini.gather;soft(:,ids)=so.gather;supports(ids)=ini.support;
 E(ids)=so.group_energy;V(ids)=so.group_noise_energy;
end
k=P.coherence;cs=coherence_scores(R0,supports,d,any(D~=0,1),k.half_width,k.dip_periods);
s=struct('R0',R0,'supports',{supports},'E',{E},'V',{V},'groups',{cs.groups},'fallback',{cs.fallback},'W',d.W);
one=struct('lambda_min',1,'lambda_max',1);
assert(norm(apply_gain(s,cs.c,one)-soft,'fro')<=1e-10*max(norm(soft,'fro'),realmin),'apply_gain does not reproduce SoftGroupRefit.');
XC=apply_gain(s,cs.c,struct('lambda_min',k.lambda_min,'lambda_max',k.lambda_max));
cc=cell2mat(reshape(cs.c,[],1));info=struct('mean_coherence',mean(cc),'groups',numel(cc));
end
