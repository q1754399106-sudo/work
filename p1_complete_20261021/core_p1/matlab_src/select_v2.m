function out=select_v2(D,models,H,lambda,sigma)
o=v2_options(); gamma=o.gamma; if nnz(H)==0, gamma=0; end
ctx=context_v2(D,models,H,lambda,gamma); target=sqrt(size(D,1)*sum(sigma.^2));
scales=o.scale_grid(:); residuals=zeros(size(scales)); candidates=cell(size(scales)); records=cell(size(scales));
for k=1:numel(scales)
 [candidates{k},records{k}]=solve_v2(ctx,scales(k),'Direct'); residuals(k)=records{k}.data_residual;
end
below=find(residuals<=target,1,'last'); above=find(residuals>=target,1,'first'); steps=0;
if ~isempty(below) && ~isempty(above) && below<above
 lo=scales(below); hi=scales(above);
 for step=1:o.refinement_max_steps
  if min(abs(residuals/target-1))<=o.residual_tolerance, break; end
  scale=sqrt(lo*hi); [S,r]=solve_v2(ctx,scale,'Direct');
  scales(end+1)=scale; residuals(end+1)=r.data_residual; candidates{end+1}=S; records{end+1}=r; steps=step;
  if r.data_residual<target, lo=scale; else, hi=scale; end
 end
end
[~,ix]=min(abs(residuals-target)); scale=scales(ix); direct=candidates{ix};
[split,sr]=solve_v2(ctx,scale,'Split'); dr=records{ix};
out=struct('direct',direct,'split',split,'scale',scale,'scale_grid',scales,'residual_grid',residuals, ...
 'target',target,'endpoint',scale==o.scale_grid(1) || scale==o.scale_grid(end), ...
 'residual_relative_mismatch',abs(dr.data_residual/target-1),'records',{records}, ...
 'direct_record',dr,'split_record',sr,'selected_index',ix,'refinement_steps',steps, ...
 'selected_method','CombinedSplit','gather',split);
end
