function [S,r]=p1_solve(ctx,scale,start)
% The covariance and pilot are fixed throughout this quadratic optimization.
if nargin<3 || isempty(start), start=zeros(ctx.n*ctx.nx,1); end
lambda=scale*ctx.ridge; mu=scale*ctx.mu;
originalOp=@(x)p1_forward(ctx.G,p1_precision(ctx.C,p1_forward(ctx.G,x,false)),true)+lambda.*x+mu*(ctx.HH*x);
if isfield(ctx,'normal_low_rank')
 op=@(x)p1_fast_normal(ctx,x,lambda,mu);
else
 op=originalOp;
end
if strcmp(ctx.options.preconditioner,'block')
 factors=cell(1,ctx.nx);
 for j=1:ctx.nx
  lamj=lambda((j-1)*ctx.n+1);
  B=ctx.block_grams{j}+lamj*eye(ctx.n)+mu*ctx.HHblocks{j};
  factors{j}=chol((B+B')/2,'lower');
 end
 pre=@(x)p1_block_precondition(x,factors,ctx.n,ctx.nx);
else
 diagPre=ctx.gram_diagonal+lambda+mu*ctx.HHdiagonal;
 pre=@(x)x./max(diagPre,realmin);
end
[s,flag,relres,iter]=pcg(op,ctx.rhs,ctx.options.pcg_tolerance,ctx.options.pcg_max_iterations,pre,[],start(:));
initialflag=flag; solver='PCG_GramWoodbury';
fastNormal=norm(op(s)-ctx.rhs)/max(norm(ctx.rhs),realmin);
% Independent original G/precision audit is performed for EVERY searched solve.
normal=norm(originalOp(s)-ctx.rhs)/max(norm(ctx.rhs),realmin);
% One stricter independent restart makes a failed warm start visible and useful.
if flag~=0 || normal>2e-7
 [s2,flag2,relres2,iter2]=pcg(originalOp,ctx.rhs,min(ctx.options.pcg_tolerance,1e-9),2*ctx.options.pcg_max_iterations,pre,[],s);
 normal2=norm(originalOp(s2)-ctx.rhs)/max(norm(ctx.rhs),realmin);
 if isfinite(normal2) && normal2<normal
  s=s2; flag=flag2; relres=relres2; iter=iter+iter2; normal=normal2; solver='PCG_OriginalPrecision_Restart';
 end
end
e=p1_forward(ctx.G,s,false)-ctx.y; pe=p1_precision(ctx.C,e); weighted=sum(e.*pe);
negativeQuadratic=weighted < -max(100*eps*norm(e)*norm(pe),realmin);
if negativeQuadratic || ~all(isfinite(s)) || ~isfinite(normal) || ~isfinite(weighted), flag=5; end
S=reshape(s,ctx.n,ctx.nx);
r=struct('flag',flag,'initial_flag',initialflag,'solver',solver,'iterations',iter, ...
 'normal_residual',normal,'data_residual',sqrt(max(weighted,0)), ...
 'pcg_relative_residual',relres,'scale',scale,'preconditioner',ctx.options.preconditioner, ...
 'mahalanobis_quadratic',weighted,'negative_precision_quadratic',negativeQuadratic, ...
 'fast_normal_residual_before_restart',fastNormal,'original_precision_audit',true);
end
function y=p1_block_precondition(x,L,n,nx)
y=zeros(size(x));
for j=1:nx, ids=(j-1)*n+(1:n); y(ids)=L{j}'\(L{j}\x(ids)); end
end
function y=p1_fast_normal(ctx,x,lambda,mu)
y=zeros(size(x));
for j=1:ctx.nx
 ids=(j-1)*ctx.n+(1:ctx.n); y(ids)=ctx.block_grams{j}*x(ids);
end
if ctx.C.rank>0
 z=ctx.normal_low_rank'*x;
 y=y-ctx.normal_low_rank*(ctx.C.T'\(ctx.C.T\z));
end
y=y+lambda.*x+mu*(ctx.HH*x);
end
