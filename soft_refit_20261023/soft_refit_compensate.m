function out=soft_refit_compensate(D,cm,offset,sigma)
% Public whole-gather entry; fixed-size batches bound dictionary memory.
% cm is the same calibrated configuration as FullCalRefit.
% Optional sigma uses the existing data-only WindowMAD estimator (0.02-0.25s).
% Users should supply their own data-based sigma when that window has signal.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'source'));
[n,nx]=size(D);assert(n==numel(cm.t)&&all(isfinite(D(:))));
if isscalar(offset),offset=repmat(offset,1,nx);end
if nargin<4
 sigma=zeros(1,nx);
 for j=1:nx
  if all(D(:,j)==0),sigma(j)=realmin;
  else,ni=estimate_trace_snr(D(:,j),cm.t);sigma(j)=ni.sigma;end
 end
end
sigma=reshape(sigma,1,[]);offset=reshape(offset,1,[]);
assert(numel(offset)==nx&&numel(sigma)==nx&&all(isfinite(sigma))&&all(sigma>0));
out=struct('gather',zeros(n,nx),'reflectivity',zeros(n,nx), ...
 'unpruned_gather',zeros(n,nx),'group_gain',{cell(1,nx)}, ...
 'optimization_flags',0,'initial_normal_residual',0,'selected_method','SoftGroupRefit');
for first=1:11:nx
 ids=first:min(first+10,nx);dicts=cell(1,numel(ids));
 for q=1:numel(ids)
  prior=find(offset(ids(1:q-1))==offset(ids(q)),1);
  if isempty(prior)
   d=prune_dictionary(cm,offset(ids(q)));Wscaled=d.W./reshape(d.scale,1,[]);
   d.source_gram=Wscaled'*Wscaled;dicts{q}=d;
  else,dicts{q}=dicts{prior};end
 end
 initial=reflectivity_refit_prepared(D(:,ids),dicts,sigma(ids));
 result=reflectivity_soft_prepared(D(:,ids),dicts,sigma(ids),initial);
 out.gather(:,ids)=result.gather;out.reflectivity(:,ids)=result.reflectivity;
 out.unpruned_gather(:,ids)=initial.gather;out.group_gain(ids)=result.group_gain;
 out.optimization_flags=out.optimization_flags+initial.optimization_flags;
 out.initial_normal_residual=max(out.initial_normal_residual,result.initial_normal_residual);
end
out.rule='gain=source group energy/(source group energy+conditional propagated white noise energy)';
out.uncertainty_assumption='Fixed selected support; white noise proxy, no Q-model bias or selection uncertainty.';
end
