function run_cv_smoke(outdir)
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_src'));
% Independent dense oracle for the coupled normal-equation solver.
old=rng;rng(20262201);n=13;nx=3;H=randn(27,n*nx);HH=H'*H;
c=struct('n',n,'nx',nx,'mu_base',.7,'HH',HH,'rhsS',randn(n,nx));
for j=1:nx
 A=randn(n);B=randn(n);c.N0S{j}=A'*A+eye(n);c.N1S{j}=B'*B+eye(n);
 ids=(j-1)*n+(1:n);c.HHblocks{j}=HH(ids,ids);
end
[S,r]=cv_solve(c,.3,.4);N=cell(1,nx);for j=1:nx,N{j}=c.N0S{j}+.3*c.N1S{j};end
exact=(blkdiag(N{:})+.4*c.mu_base*HH)\c.rhsS(:);relative=norm(S(:)-exact)/norm(exact);assert(relative<2e-6&&r.flag==0);
% Check that masking occurs before nonlocal whitening.
D=randn(601,11);masks=cv_masks(601,11);K=randn(601);mask=masks{1};perturb=randn(size(D)).*~mask;
leak=norm(K*(D.*mask)-K*((D+perturb).*mask),'fro');assert(leak==0);
rng(old);run_cv_screen(outdir,1);
uq_writejson(fullfile(outdir,'smoke_checks.json'),struct('dense_oracle_relative_error',relative,'mask_before_K_leakage',leak,'passed',true));
end
