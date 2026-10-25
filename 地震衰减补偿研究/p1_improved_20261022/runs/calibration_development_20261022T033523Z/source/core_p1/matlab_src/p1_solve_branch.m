function [out,prior,preset]=p1_solve_branch(D,models,map,cm,domain,R_avo,branch_name)
% Unified P1 branch solver for Direct and Split routes with identical prior and error budget.
root=fileparts(fileparts(mfilename('fullpath')));
p1spec=jsondecode(fileread(fullfile(root,'inputs','p1_experiment_spec.json')));
p1=p1spec.p1_hyperparameters;
o_struct=structure_options(); v2=v2_options();
[n,nx]=size(D); sigma=zeros(1,nx); lam=sigma; pilot=zeros(n,nx);
for j=1:nx
 ni=estimate_trace_snr(D(:,j),cm.t); sigma(j)=ni.sigma;
 [~,lam(j)]=discrepancy_apply(models{j}.direct,D(:,j),cm.lambda_grid,sqrt(n)*sigma(j));
 F=models{j}.direct; u=F.U'*D(:,j); pilot(:,j)=F.V*(F.z./(F.z.^2+4*lam(j)).*u);
end
if any(strcmp(branch_name,{'Baseline','Frozen','StrongerDamping'}))
 gamma=10; if strcmp(branch_name,'Baseline'), gamma=0; end
 eta=1.0; if strcmp(branch_name,'StrongerDamping'), eta=1.15; end
 preset=struct('name',branch_name,'gamma',gamma,'pilot_multiplier',4,'noise_target_multiplier',eta);
 prior=uq_build_prior(D,models,cm,domain,preset);
 out=uq_select(D,models,prior,preset);
 return;
end
use_noise=any(strcmp(branch_name,{'NoiseUncertaintyOnly','CombinedUncertainty','CombinedUncertainty_NoStructure'}));
use_q=any(strcmp(branch_name,{'QModelErrorOnly','CombinedUncertainty','CombinedUncertainty_NoStructure'}));
use_struct=~strcmp(branch_name,'CombinedUncertainty_NoStructure');
if use_noise
 [D_work,b_coh,ortho_ratio,tail_mask,extra_diag]=p1_noise_uncertainty(D,pilot,models,sigma,lam,domain,R_avo,p1);
else
 D_work=D; b_coh=zeros(n,nx); ortho_ratio=0;
 [~,~,~,tail_mask,extra_diag]=p1_noise_uncertainty(D,pilot,models,sigma,lam,domain,R_avo,p1);
end
if ~use_struct
 tail_mask=false(n,nx); extra_diag=zeros(n,nx);
end
pilot_w=pilot;
if any(b_coh(:)~=0)
 for j=1:nx
  F=models{j}.direct; u=F.U'*D_work(:,j); pilot_w(:,j)=F.V*(F.z./(F.z.^2+4*lam(j)).*u);
 end
end
gamma=10*double(use_struct);
if use_struct
 [H0,~]=learn_structure_operator(pilot_w,cm,o_struct);
 if strcmp(domain,'Offset')
  w_avo=p1.offset_avo_ortho_weight_qonly; if use_noise, w_avo=p1.offset_avo_ortho_weight_noise; end
  H=[H0*R_avo; w_avo*R_avo];
 else
  if any(tail_mask(:))
   row_scale=ones(n,nx); row_scale(tail_mask)=p1.section_tail_h_scale;
   H=spdiags(row_scale(:),0,n*nx,n*nx)*H0;
  else
   H=H0;
  end
 end
else
 H=sparse(n*nx,n*nx);
end
[r2_D,r2_S,eta_eff,b_eD,b_eS]=p1_q_model_error(pilot_w,models,cm,domain,use_q,use_noise,p1);
D_D=D_work-b_eD; D_S=D_work-b_eS;
Yw=zeros(n,nx); rhsD=Yw; rhsS=Yw; uD=Yw; uS=Yw;
for j=1:nx
 m=models{j}; Yw(:,j)=m.K*D_S(:,j);
 rhsD(:,j)=m.G'*D_D(:,j); rhsS(:,j)=m.F'*Yw(:,j);
 uD(:,j)=m.direct.U'*D_D(:,j); uS(:,j)=m.split.U'*Yw(:,j);
end
HH=H'*H; blocks=cell(1,nx);
for j=1:nx, ids=(j-1)*n+(1:n); blocks{j}=full(HH(ids,ids)); end
ctx=struct('D',D_D,'Yw',Yw,'models',{models},'H',H,'HH',HH,'HHblocks',{blocks}, ...
 'lambda',lam,'mu_base',gamma*median(lam),'rhsD',rhsD,'rhsS',rhsS,'uD',uD,'uS',uS, ...
 'r2_D',{r2_D},'r2_S',{r2_S},'extra_diag',extra_diag,'n',n,'nx',nx);
% Exact affine pencil: N(scale)=N0+scale*N1. Cache scale-invariant dense
% products once per observation/branch instead of rebuilding at every scale.
% Same normal equations; floating-point regrouping is validated by replay.
ctx.N0D=cell(1,nx); ctx.N0S=cell(1,nx); ctx.N1D=cell(1,nx); ctx.N1S=cell(1,nx); ctx.exact_final=false;
for j=1:nx
 for ir=1:2
  if ir==1, F=models{j}.direct; r2=r2_D{j}; else, F=models{j}.split; r2=r2_S{j}; end
  N0=(F.V.*reshape(F.z.^2,1,[]))*F.V';
  N1=(F.V.*reshape(lam(j)*r2,1,[]))*F.V'+diag(extra_diag(:,j));
  if ir==1, ctx.N0D{j}=N0; ctx.N1D{j}=N1; else, ctx.N0S{j}=N0; ctx.N1S{j}=N1; end
 end
end
target=eta_eff*sqrt(n*sum(sigma.^2));
scales=v2.scale_grid(:); residuals=zeros(size(scales)); candidates=cell(size(scales)); records=cell(size(scales));
for k=1:numel(scales)
 [candidates{k},records{k}]=solve_weighted(ctx,scales(k),'Direct',o_struct); residuals(k)=records{k}.data_residual;
end
below=find(residuals<=target,1,'last'); above=find(residuals>=target,1,'first');
if ~isempty(below) && ~isempty(above) && below<above
 lo=scales(below); hi=scales(above);
 for step=1:v2.refinement_max_steps
  if min(abs(residuals/target-1))<=v2.residual_tolerance, break; end
  sc=sqrt(lo*hi); [S,r]=solve_weighted(ctx,sc,'Direct',o_struct);
  scales(end+1)=sc; residuals(end+1)=r.data_residual; candidates{end+1}=S; records{end+1}=r;
  if r.data_residual<target, lo=sc; else, hi=sc; end
 end
end
[~,ix]=min(abs(residuals-target)); scale=scales(ix);
% Rebuild the selected normal equations in the original arithmetic order.
% The cache speeds scale selection; reported reconstructions use original N.
ctx.exact_final=true; [direct,dr]=solve_weighted(ctx,scale,'Direct',o_struct); records{ix}=dr;
[split,sr]=solve_weighted(ctx,scale,'Split',o_struct);
prior=struct('H',H,'sigma',sigma,'lambda',lam,'pilot',pilot_w,'ortho_ratio',ortho_ratio);
preset=struct('name',branch_name,'gamma',gamma,'pilot_multiplier',4,'noise_target_multiplier',eta_eff);
out=struct('direct',direct,'split',split,'scale',scale,'scale_grid',scales,'residual_grid',residuals,'target',target, ...
 'endpoint',scale==v2.scale_grid(1) || scale==v2.scale_grid(end),'residual_mismatch',abs(residuals(ix)/target-1), ...
 'records',{[records(:);{sr}]},'direct_record',records{ix},'split_record',sr);
end
function [S,r]=solve_weighted(ctx,scale,route,o)
n=ctx.n; nx=ctx.nx; lam=scale*ctx.lambda; mu=scale*ctx.mu_base;
N=cell(1,nx); L=N; S=zeros(n,nx);
if strcmp(route,'Direct'), rhs=ctx.rhsD; u=ctx.uD; r2_list=ctx.r2_D; else, rhs=ctx.rhsS; u=ctx.uS; r2_list=ctx.r2_S; end
for j=1:nx
 m=ctx.models{j}; if strcmp(route,'Direct'), F=m.direct; else, F=m.split; end
 r2=r2_list{j}; z=F.z; z2=z.^2;
 d_mode=z2+lam(j)*r2;
 if ctx.exact_final
  Nj=F.V*diag(d_mode)*F.V'+diag(scale*ctx.extra_diag(:,j));
 elseif strcmp(route,'Direct'), Nj=ctx.N0D{j}+scale*ctx.N1D{j};
 else, Nj=ctx.N0S{j}+scale*ctx.N1S{j}; end
 N{j}=Nj;
 if any(ctx.extra_diag(:,j)~=0) && nnz(ctx.H)==0
  Lj=chol(Nj,'lower'); S(:,j)=Lj'\(Lj\rhs(:,j));
 else
  S(:,j)=F.V*(z./d_mode.*u(:,j));
 end
 if nnz(ctx.H)>0, L{j}=chol(Nj+mu*ctx.HHblocks{j},'lower'); end
end
flag=0; iter=0; initialflag=0; solver='SVD';
if nnz(ctx.H)>0
 op=@(v)prod_b(v,N,n,nx)+mu*(ctx.HH*v); pre=@(v)sol_b(v,L,n,nx);
 [s,flag,~,iter]=pcg(op,rhs(:),o.pcg_tolerance,o.pcg_max_iterations,pre,[],S(:)); initialflag=flag; solver='PCG';
 if flag~=0 || norm(op(s)-rhs(:))/norm(rhs(:))>2*o.pcg_tolerance
  cells=cellfun(@sparse,N,'UniformOutput',false); B=blkdiag(cells{:})+mu*ctx.HH;
  [C,p]=chol(B,'lower'); assert(p==0); s=C'\(C\rhs(:)); flag=0; solver='SparseCholesky';
 end
 normal=norm(op(s)-rhs(:))/norm(rhs(:)); S=reshape(s,n,nx);
else
 check=zeros(n,nx); for j=1:nx, check(:,j)=N{j}*S(:,j)-rhs(:,j); end
 normal=norm(check(:))/norm(rhs(:));
end
rr=0;
for j=1:nx
 if strcmp(route,'Direct'), e=ctx.models{j}.G*S(:,j)-ctx.D(:,j); else, e=ctx.models{j}.F*S(:,j)-ctx.Yw(:,j); end
 rr=rr+sum(e.^2);
end
r=struct('flag',flag,'initial_flag',initialflag,'solver',solver,'iterations',iter,'normal_residual',normal,'data_residual',sqrt(rr));
end
function y=prod_b(v,N,n,nx), v=reshape(v,n,nx); y=zeros(n,nx); for j=1:nx, y(:,j)=N{j}*v(:,j); end; y=y(:); end
function y=sol_b(v,L,n,nx), v=reshape(v,n,nx); y=zeros(n,nx); for j=1:nx, y(:,j)=L{j}'\(L{j}\v(:,j)); end; y=y(:); end
