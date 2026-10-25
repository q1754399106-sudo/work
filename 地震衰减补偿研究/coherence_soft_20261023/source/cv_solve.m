function [S,record]=cv_solve(c,alpha,beta_factor)
% Split estimator with independent alpha and beta, fixed per solve.
o=structure_options();n=c.n;nx=c.nx;beta=beta_factor*c.mu_base;
N=cell(1,nx);L=N;
for j=1:nx
 N{j}=c.N0S{j}+alpha*c.N1S{j};N{j}=(N{j}+N{j}')/2;
 L{j}=chol(N{j}+beta*c.HHblocks{j},'lower');
end
op=@(x)prodblocks(x,N,n,nx)+beta*(c.HH*x);
pre=@(x)solveblocks(x,L,n,nx);
[s,flag,~,iter]=pcg(op,c.rhsS(:),o.pcg_tolerance,o.pcg_max_iterations,pre);
initial=flag;normal=norm(op(s)-c.rhsS(:))/max(norm(c.rhsS(:)),eps);solver='PCG';
if flag~=0||normal>2*o.pcg_tolerance
 B=blkdiag(N{:})+beta*c.HH;C=chol(B,'lower');s=C'\(C\c.rhsS(:));flag=0;solver='Cholesky';
 normal=norm(op(s)-c.rhsS(:))/max(norm(c.rhsS(:)),eps);
end
S=reshape(s,n,nx);record=struct('flag',flag,'initial_flag',initial,'normal_residual',normal,'iterations',iter,'solver',solver,'alpha',alpha,'beta_factor',beta_factor,'beta',beta);
end
function y=prodblocks(x,N,n,nx)
x=reshape(x,n,nx);y=zeros(n,nx);for j=1:nx,y(:,j)=N{j}*x(:,j);end;y=y(:);
end
function y=solveblocks(x,L,n,nx)
x=reshape(x,n,nx);y=zeros(n,nx);for j=1:nx,y(:,j)=L{j}'\(L{j}\x(:,j));end;y=y(:);
end
