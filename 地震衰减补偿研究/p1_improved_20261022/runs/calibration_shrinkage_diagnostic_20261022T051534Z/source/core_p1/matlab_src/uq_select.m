function out=uq_select(D,models,prior,preset)
% Same noise-target selection as v2, explicit preset overrides, no true Q.
o=v2_options(); gamma=preset.gamma; if nnz(prior.H)==0, gamma=0; end
ctx=context_v2(D,models,prior.H,prior.lambda,gamma);
target=preset.noise_target_multiplier*sqrt(size(D,1)*sum(prior.sigma.^2));
scales=o.scale_grid(:); residuals=zeros(size(scales)); candidates=cell(size(scales)); records=cell(size(scales));
for k=1:numel(scales), [candidates{k},records{k}]=solve_v2(ctx,scales(k),'Direct'); residuals(k)=records{k}.data_residual; end
below=find(residuals<=target,1,'last'); above=find(residuals>=target,1,'first');
if ~isempty(below) && ~isempty(above) && below<above
 lo=scales(below); hi=scales(above);
 for step=1:o.refinement_max_steps
  if min(abs(residuals/target-1))<=o.residual_tolerance, break; end
  scale=sqrt(lo*hi); [S,r]=solve_v2(ctx,scale,'Direct');
  scales(end+1)=scale; residuals(end+1)=r.data_residual; candidates{end+1}=S; records{end+1}=r;
  if r.data_residual<target, lo=scale; else, hi=scale; end
 end
end
[~,ix]=min(abs(residuals-target)); scale=scales(ix); direct=candidates{ix}; [split,sr]=solve_v2(ctx,scale,'Split');
out=struct('direct',direct,'split',split,'scale',scale,'scale_grid',scales,'residual_grid',residuals,'target',target, ...
 'endpoint',scale==o.scale_grid(1) || scale==o.scale_grid(end),'residual_mismatch',abs(residuals(ix)/target-1), ...
 'records',{[records(:);{sr}]},'direct_record',records{ix},'split_record',sr);
end
