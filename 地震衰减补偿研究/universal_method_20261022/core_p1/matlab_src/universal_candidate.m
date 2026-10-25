function out=universal_candidate(ctx,cm,profile)
% Fixed rules for every geology and acquisition domain. No truth/scene input.
assert(any(strcmp(profile,{'SourceMixture','SourceOnly'})));
P=p1_source_precision(ctx.n,cm.dt,cm.fpeak,.02);
blend=.5;
if strcmp(profile,'SourceOnly'),blend=1;end
for j=1:ctx.nx
 E=diag(ctx.extra_diag(:,j));
 ctx.N1D{j}=(1-blend)*(ctx.N1D{j}-E)+blend*ctx.lambda(j)*P+E;
 ctx.N1S{j}=(1-blend)*(ctx.N1S{j}-E)+blend*ctx.lambda(j)*P+E;
end
% SourceOnly also ablates the learned structure term, globally.
beta=.3;if strcmp(profile,'SourceOnly'),beta=0;end
out=joint_select(ctx,'DecoupledWeak',beta,{});
out.profile=profile;out.source_blend=blend;out.source_floor=.02;
end
