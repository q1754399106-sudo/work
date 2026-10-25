function [g,diag]=run_scene(cm,D,sigma,referenceQ,errorpct,opts)
% All methods share observation, sigma, wavelet and the applied (possibly wrong) Q.
% Truth is never an input. Blocks of 11 traces as in FullCalRefit.
offset=450;if isfield(cm,'solver_offset'),offset=cm.solver_offset;end   % independent scenes: 0
used=cm;used.Q=reshape(referenceQ,1,[])*(1+errorpct/100);d=prune_dictionary(used,offset);
Ws=d.W./reshape(d.scale,1,[]);d.source_gram=Ws'*Ws;nx=size(D,2);sigma=reshape(sigma,1,[]);
lean=isfield(opts,'lean')&&opts.lean;   % frozen runs: FullCalRefit, Soft and QEnsembleMean only
z=zeros(size(D));g=struct('FullCalRefit',z,'SoftGroupRefit',z);if ~lean,g.PrunedFullCalRefit=z;end
for s=opts.s_grid
 g.(tag('QEnsembleMean',s))=z;if ~lean,g.(tag('BAE',s))=z;g.(tag('BAECovOnly',s))=z;end
end
ops=struct();ens=struct();
for s=opts.s_grid
 if ~lean,ops.(tag('s',s))=q_sample_operators(used,offset,s,opts.M,opts.rho,opts.seed);end
 nodes=sqrt(2)*erfinv(2*((1:opts.ensemble_nodes)-.5)/opts.ensemble_nodes-1);dd=cell(1,numel(nodes));
 for k=1:numel(nodes),c=used;c.Q=used.Q*exp(s*nodes(k));dd{k}=prune_dictionary(c,offset);end
 ens.(tag('s',s))=dd;
end
flags=0;normal=0;rank=zeros(1,nx);
for first=1:11:nx
 ids=first:min(first+10,nx);dicts=repmat({d},1,numel(ids));
 initial=reflectivity_refit_prepared(D(:,ids),dicts,sigma(ids));
 soft=reflectivity_soft_prepared(D(:,ids),dicts,sigma(ids),initial);
 g.FullCalRefit(:,ids)=initial.gather;g.SoftGroupRefit(:,ids)=soft.gather;
 flags=flags+initial.optimization_flags;normal=max(normal,initial.normal_residual);
 if ~lean
  pruned=reflectivity_prune_prepared(D(:,ids),dicts,sigma(ids),initial);g.PrunedFullCalRefit(:,ids)=pruned.gather;
  flags=flags+pruned.optimization_flags;normal=max(normal,pruned.normal_residual);
 end
 for s=opts.s_grid
  if ~lean
   Gs=ops.(tag('s',s));
   b=bae_q_refit(D(:,ids),d,Gs,sigma(ids),initial,true);g.(tag('BAE',s))(:,ids)=b.gather;rank(ids)=b.error_rank;
   c=bae_q_refit(D(:,ids),d,Gs,sigma(ids),initial,false);g.(tag('BAECovOnly',s))(:,ids)=c.gather;
   flags=flags+b.optimization_flags+c.optimization_flags;normal=max([normal,b.normal_residual,c.normal_residual]);
  end
  dd=ens.(tag('s',s));acc=zeros(size(D,1),numel(ids));
  for k=1:numel(dd)
   e=reflectivity_refit_prepared(D(:,ids),repmat(dd(k),1,numel(ids)),sigma(ids));acc=acc+e.gather/numel(dd);
   flags=flags+e.optimization_flags;normal=max(normal,e.normal_residual);
  end
  g.(tag('QEnsembleMean',s))(:,ids)=acc;
 end
end
diag=struct('used_Q',used.Q,'optimization_flags',flags,'max_normal_residual',normal,'max_error_rank',max(rank));
end
function name=tag(prefix,s)
name=sprintf('%s_s%03d',prefix,round(100*s));
end
