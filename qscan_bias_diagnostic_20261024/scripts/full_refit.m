function [S,info]=full_refit(D,sigma,cm,Q,offset)
% FullCalRefit with applied Q, blocks of 11 traces as in run_scene.m and qscan_refit.m.
c=cm;c.Q=Q;d=prune_dictionary(c,offset);[n,nx]=size(D);S=zeros(n,nx);flags=0;normal=0;support=zeros(1,nx);
for first=1:11:nx
 ids=first:min(first+10,nx);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sigma(ids));
 S(:,ids)=o.gather;flags=flags+o.optimization_flags;normal=max(normal,o.normal_residual);support(ids)=o.support_after;
end
info=struct('optimization_flags',flags,'max_normal_residual',normal,'support_size',support);
end
