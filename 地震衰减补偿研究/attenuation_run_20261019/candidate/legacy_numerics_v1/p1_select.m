function out=p1_select(directCtx,splitCtx,pilot)
% Select only on the Direct whitened discrepancy; matched scale for Split.
o=directCtx.options; target=sqrt(directCtx.C.target_squared);
scales=o.scale_grid(:); residuals=zeros(size(scales)); candidates=cell(size(scales)); records=cell(size(scales));
start=pilot;
for k=1:numel(scales)
 [candidates{k},records{k}]=p1_solve(directCtx,scales(k),start);
 residuals(k)=records{k}.data_residual; start=candidates{k};
end
below=find(residuals<=target,1,'last'); above=find(residuals>=target,1,'first'); steps=0;
if ~isempty(below) && ~isempty(above) && below<above
 lo=scales(below); hi=scales(above);
 for step=1:o.refinement_max_steps
  if min(abs(residuals/target-1))<=o.residual_tolerance, break; end
  scale=sqrt(lo*hi); [S,r]=p1_solve(directCtx,scale,start); start=S;
  scales(end+1)=scale; residuals(end+1)=r.data_residual; candidates{end+1}=S; records{end+1}=r; steps=step;
  if r.data_residual<target, lo=scale; else, hi=scale; end
 end
end
[~,ix]=min(abs(residuals-target)); scale=scales(ix); direct=candidates{ix};
[split,sr]=p1_solve(splitCtx,scale,direct);
out=struct('direct',direct,'split',split,'scale',scale,'scale_grid',scales,'residual_grid',residuals, ...
 'target',target,'split_target',sqrt(splitCtx.C.target_squared), ...
 'endpoint',scale==o.scale_grid(1) || scale==o.scale_grid(end), ...
 'residual_mismatch',abs(residuals(ix)/target-1),'records',{[records(:);{sr}]}, ...
 'direct_record',records{ix},'split_record',sr,'selected_index',ix,'refinement_steps',steps, ...
 'residual_units','fixed-covariance Mahalanobis norm', ...
 'residual_monotonicity_violations',sum(diff(residuals(1:numel(o.scale_grid))) < -1e-7*max(target,1)));
end
