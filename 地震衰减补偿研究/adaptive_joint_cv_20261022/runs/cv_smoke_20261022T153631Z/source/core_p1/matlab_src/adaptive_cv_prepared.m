function out=adaptive_cv_prepared(D,models,map,cm,domain,R,training_models)
% One global rule. No model/scenario labels, clean traces or true Q inputs.
root=fileparts(fileparts(mfilename('fullpath')));
plan=jsondecode(fileread(fullfile(root,'inputs','adaptive_cv_plan.json')));
[n,nx]=size(D);masks=cv_masks(n,nx);alphas=plan.alpha_grid(:);betas=plan.beta_grid(:);
loss=zeros(numel(alphas),numel(betas),numel(masks));normal=0;flags=0;
cm.joint_prepare_only=true;
for fold=1:numel(masks)
 mask=masks{fold};trainD=D.*mask;traincm=cm;traincm.training_mask=mask;
 if nargin<7||isempty(training_models),tm=cv_training_models(models,mask);else,tm=training_models{fold};end
 [~,prior,~,ctx]=p1_solve_branch(trainD,tm,map,traincm,domain,R,'CombinedUncertainty');ctx=cv_normalize_context(ctx);
 for ia=1:numel(alphas)
  for ib=1:numel(betas)
   [S,r]=cv_solve(ctx,alphas(ia),betas(ib));normal=max(normal,r.normal_residual);flags=flags+r.flag;
   value=0;count=0;
   for j=1:nx
    held=~mask(:,j);e=models{j}.G*S(:,j)-D(:,j);
    value=value+sum(e(held).^2)/max(prior.sigma(j)^2,realmin);count=count+sum(held);
   end
   loss(ia,ib,fold)=value/count;
  end
 end
end
average=mean(loss,3);[~,ix]=min(average(:));[ia,ib]=ind2sub(size(average),ix);
[~,prior,preset,ctx]=p1_solve_branch(D,models,map,cm,domain,R,'CombinedUncertainty');ctx=cv_normalize_context(ctx);
[S,r]=cv_solve(ctx,alphas(ia),betas(ib));normal=max(normal,r.normal_residual);flags=flags+r.flag;
assert(flags==0&&normal<2e-7);
out=struct('gather',S,'alpha',alphas(ia),'beta_factor',betas(ib),'beta',r.beta,'cv_loss',loss,'cv_mean',average, ...
 'alpha_grid',alphas,'beta_grid',betas,'selection_index',[ia,ib],'alpha_endpoint',ia==1||ia==numel(alphas), ...
 'beta_endpoint',ib==1||ib==numel(betas),'final_record',r,'max_normal',normal,'flags',flags,'prior',prior,'preset',preset, ...
 'selected_method','AdaptiveJointCV_Split','development_only',true);
end
