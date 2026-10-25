function [S,r]=solve_v2(ctx,scale,route)
o=structure_options(); n=ctx.n; nx=ctx.nx; lam=scale*ctx.lambda; mu=scale*ctx.mu_base;
N=cell(1,nx); L=N; S=zeros(n,nx);
if strcmp(route,'Direct'), rhs=ctx.rhsD; u=ctx.uD; else, rhs=ctx.rhsS; u=ctx.uS; end
for j=1:nx
 m=ctx.models{j};
 if strcmp(route,'Direct'), F=m.direct; gram=m.gram_direct; else, F=m.split; gram=m.gram_split; end
 S(:,j)=F.V*(F.z./(F.z.^2+lam(j)).*u(:,j));
 if nnz(ctx.H)>0
  N{j}=gram+lam(j)*eye(n); L{j}=chol(N{j}+mu*ctx.HHblocks{j},'lower');
 end
end
flag=0; iter=0; initialflag=0; solver='SVD';
if nnz(ctx.H)>0
 op=@(v)product(v,N,n,nx)+mu*(ctx.HH*v); pre=@(v)solve(v,L,n,nx);
 [s,flag,~,iter]=pcg(op,rhs(:),o.pcg_tolerance,o.pcg_max_iterations,pre,[],S(:)); initialflag=flag; solver='PCG';
 if flag~=0 || norm(op(s)-rhs(:))/norm(rhs(:))>2*o.pcg_tolerance
  cells=cellfun(@sparse,N,'UniformOutput',false); B=blkdiag(cells{:})+mu*ctx.HH;
  [C,p]=chol(B,'lower'); assert(p==0); s=C'\(C\rhs(:)); flag=0; solver='SparseCholesky';
 end
 normal=norm(op(s)-rhs(:))/norm(rhs(:)); S=reshape(s,n,nx);
else
 check=zeros(n,nx);
 for j=1:nx
  m=ctx.models{j}; if strcmp(route,'Direct'), gram=m.gram_direct; else, gram=m.gram_split; end
  check(:,j)=gram*S(:,j)+lam(j)*S(:,j)-rhs(:,j);
 end
 normal=norm(check(:))/norm(rhs(:));
end
rr=0;
for j=1:nx
 if strcmp(route,'Direct'), e=ctx.models{j}.G*S(:,j)-ctx.D(:,j);
 else, e=ctx.models{j}.F*S(:,j)-ctx.Yw(:,j); end
 rr=rr+sum(e.^2);
end
r=struct('flag',flag,'initial_flag',initialflag,'solver',solver,'iterations',iter,'normal_residual',normal,'data_residual',sqrt(rr));
end
function y=product(v,N,n,nx)
v=reshape(v,n,nx); y=zeros(n,nx); for j=1:nx, y(:,j)=N{j}*v(:,j); end; y=y(:);
end
function y=solve(v,L,n,nx)
v=reshape(v,n,nx); y=zeros(n,nx); for j=1:nx, y(:,j)=L{j}'\(L{j}\v(:,j)); end; y=y(:);
end
