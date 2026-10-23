function test_joint_algebra(outdir)
% Independent dense normal-equation oracle; uses synthetic small matrices.
root=fileparts(mfilename('fullpath'));oldpath=path;oldrng=rng;cleanup=onCleanup(@()restore(oldpath,oldrng));
addpath(fullfile(root,'matlab_src'));assert(~isfolder(outdir));mkdir(outdir);rng(701);
n=8;nx=3;D=randn(n,nx);H=randn(7,n*nx);HH=H'*H;models=cell(1,nx);ensemble=models;
c=struct('n',n,'nx',nx,'H',H,'HH',HH,'HHblocks',{cell(1,nx)},'lambda',[.1 .2 .3], ...
 'mu_base',.4,'D',D,'Yw',D,'models',{models},'target',1,'extra_diag',.01*ones(n,nx));
for j=1:nx
 G=randn(n)/sqrt(n);K=diag(linspace(.7,1,n));F=K*G;B=randn(n)/10;
 c.models{j}=struct('G',G,'K',K,'F',F);c.Yw(:,j)=K*D(:,j);
 c.N0D{j}=G'*G;c.N0S{j}=F'*F;
 c.N1D{j}=c.lambda(j)*eye(n)+diag(c.extra_diag(:,j));c.N1S{j}=c.N1D{j};
 c.rhsD(:,j)=G'*D(:,j);c.rhsS(:,j)=F'*c.Yw(:,j);
 ids=(j-1)*n+(1:n);c.HHblocks{j}=HH(ids,ids);
 ensemble{j}=struct('Gbar',G,'RQ',B'*B);
end
maxerr=0;hscale=norm(H,'fro')/sqrt(size(H,1));
for method={'Decoupled','SharedIsotropic','Physical','CoupledReplay'}
 mode=method{1};out=joint_select(c,mode,.3,ensemble);
 for route={'Direct','Split'}
  name=route{1};N=cell(1,nx);rhs=c.rhsD;actual=out.direct;
  if strcmp(name,'Split'),rhs=c.rhsS;actual=out.split;end
  for j=1:nx
   if strcmp(name,'Direct'),N{j}=c.N0D{j}+out.scale*c.N1D{j};else,N{j}=c.N0S{j}+out.scale*c.N1S{j};end
   if strcmp(mode,'Physical'),N{j}=N{j}+ensemble{j}.RQ;end
  end
  beta=.3*c.mu_base/hscale^2;if strcmp(mode,'CoupledReplay'),beta=out.scale*c.mu_base;end
  oracle=(blkdiag(N{:})+beta*HH)\rhs(:);
  err=norm(actual(:)-oracle)/norm(oracle);maxerr=max(maxerr,err);assert(err<1e-6);
 end
end
uq_writejson(fullfile(outdir,'checks.json'),struct('completed',true,'max_dense_oracle_relative_error',maxerr,'routes_tested',2,'methods_tested',4));
fprintf('Independent dense oracle max error %.3e\n',maxerr);
end
function restore(p,r),path(p);rng(r);end
