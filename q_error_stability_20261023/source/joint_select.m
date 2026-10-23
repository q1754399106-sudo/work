function out=joint_select(ctx,mode,beta_factor,ensemble)
% alpha controls modal/local damping, beta is fixed during discrepancy search.
% Physical and SharedIsotropic share lambda*I+E, H, beta and mean operator;
% only Physical adds the unscaled expectation penalty RQ with coefficient 1.
o=structure_options();v=v2_options();n=ctx.n;nx=ctx.nx;
active=sum(full(sum(ctx.H.^2,2))>0);hscale=norm(ctx.H,'fro')/sqrt(max(active,1));
if hscale==0,hscale=1;end
ctx.H=ctx.H/hscale;ctx.HH=ctx.HH/hscale^2;
for j=1:nx,ctx.HHblocks{j}=ctx.HHblocks{j}/hscale^2;end
% beta in normalized coordinates: explicit, independent of alpha.
ctx.beta=beta_factor*ctx.mu_base;
ctx.coupled_replay=strcmp(mode,'CoupledReplay');
if ctx.coupled_replay,ctx.beta=ctx.mu_base*hscale^2;end
ctx.shared=any(strcmp(mode,{'Physical','SharedIsotropic'}));
ctx.qterm=strcmp(mode,'Physical');ctx.ensemble=ensemble;
if ctx.shared
 for j=1:nx
  ctx.N1D{j}=ctx.lambda(j)*eye(n)+diag(ctx.extra_diag(:,j));ctx.N1S{j}=ctx.N1D{j};
  G=ensemble{j}.Gbar;F=ctx.models{j}.K*G;
  ctx.models{j}.G=G;ctx.models{j}.F=F;
  ctx.N0D{j}=G'*G;ctx.N0S{j}=F'*F;
  ctx.rhsD(:,j)=G'*ctx.D(:,j);ctx.rhsS(:,j)=F'*ctx.Yw(:,j);
 end
end
scales=v.scale_grid(:);res=zeros(size(scales));records=cell(size(scales));
for k=1:numel(scales),[~,records{k}]=solve(ctx,scales(k),'Direct',o);res(k)=records{k}.data_residual;end
% Search every sign-changing interval; do not assume monotonicity for
% noncommuting independent penalties. Preserve endpoint diagnostics.
for step=1:v.refinement_max_steps
 if min(abs(res/ctx.target-1))<=v.residual_tolerance,break;end
 [ss,ix]=sort(scales);rr=res(ix)-ctx.target;cross=find(rr(1:end-1).*rr(2:end)<=0);
 if isempty(cross),break;end
 [~,best]=min(abs(rr(cross))+abs(rr(cross+1)));a=cross(best);s=sqrt(ss(a)*ss(a+1));
 [~,r]=solve(ctx,s,'Direct',o);scales(end+1)=s;res(end+1)=r.data_residual;records{end+1}=r;
end
[~,ix]=min(abs(res-ctx.target));scale=scales(ix);
[direct,dr]=solve(ctx,scale,'Direct',o);[split,sr]=solve(ctx,scale,'Split',o);
regularizer_difference=0;
for j=1:nx
 regularizer_difference=max(regularizer_difference,norm(ctx.N1D{j}-ctx.N1S{j},'fro')/max(norm(ctx.N1D{j},'fro'),eps));
end
out=struct('direct',direct,'split',split,'scale',scale,'scale_grid',scales,'residual_grid',res, ...
 'target',ctx.target,'endpoint',scale==v.scale_grid(1)||scale==v.scale_grid(end), ...
 'residual_mismatch',abs(dr.data_residual/ctx.target-1),'records',{[records(:);{dr};{sr}]}, ...
 'direct_record',dr,'split_record',sr,'beta',ctx.beta,'H_rms',hscale,'shared_regularizer',ctx.shared, ...
 'regularizer_difference',regularizer_difference);
end
function [S,r]=solve(c,alpha,route,o)
n=c.n;nx=c.nx;N=cell(1,nx);L=N;rhs=c.rhsD;
if c.coupled_replay,c.beta=alpha*c.beta;end
if strcmp(route,'Split'),rhs=c.rhsS;end
for j=1:nx
 if strcmp(route,'Direct'),B=c.N0D{j};R=c.N1D{j};else,B=c.N0S{j};R=c.N1S{j};end
 N{j}=B+alpha*R;if c.qterm,N{j}=N{j}+c.ensemble{j}.RQ;end
 N{j}=(N{j}+N{j}')/2;L{j}=chol(N{j}+c.beta*c.HHblocks{j},'lower');
end
op=@(x)blockprod(x,N,n,nx)+c.beta*(c.HH*x);pre=@(x)blocksolve(x,L,n,nx);
[s,flag,~,iter]=pcg(op,rhs(:),o.pcg_tolerance,o.pcg_max_iterations,pre);
initial=flag;solver='PCG';normal=norm(op(s)-rhs(:))/max(norm(rhs(:)),eps);
if flag~=0||normal>2*o.pcg_tolerance
 B=blkdiag(N{:})+c.beta*c.HH;C=chol(B,'lower');s=C'\(C\rhs(:));flag=0;solver='Cholesky';
 normal=norm(op(s)-rhs(:))/max(norm(rhs(:)),eps);
end
S=reshape(s,n,nx);rr=0;qpen=0;modal=0;
for j=1:nx
 if strcmp(route,'Direct'),e=c.models{j}.G*S(:,j)-c.D(:,j);R=c.N1D{j};
 else,e=c.models{j}.F*S(:,j)-c.Yw(:,j);R=c.N1S{j};end
 rr=rr+sum(e.^2);modal=modal+alpha*S(:,j)'*R*S(:,j);
 if c.qterm,qpen=qpen+S(:,j)'*c.ensemble{j}.RQ*S(:,j);end
end
r=struct('flag',flag,'initial_flag',initial,'solver',solver,'iterations',iter, ...
 'normal_residual',normal,'data_residual',sqrt(rr),'modal_penalty',modal, ...
 'structure_penalty',c.beta*norm(c.H*S(:))^2,'q_penalty',qpen);
end
function y=blockprod(x,N,n,nx)
x=reshape(x,n,nx);y=zeros(n,nx);for j=1:nx,y(:,j)=N{j}*x(:,j);end;y=y(:);
end
function y=blocksolve(x,L,n,nx)
x=reshape(x,n,nx);y=zeros(n,nx);for j=1:nx,y(:,j)=L{j}'\(L{j}\x(:,j));end;y=y(:);
end
