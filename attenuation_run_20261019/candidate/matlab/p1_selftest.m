function report=p1_selftest(outdir)
% Meaningful small-matrix checks, independent of the large geological suite.
before=rng; stream=RandStream('mt19937ar','Seed',17419);
n=12; nx=2; N=n*nx; rank=3; U=.08*randn(stream,N,rank);
blocks=cell(1,nx); G=blocks; M=blocks; fullG=zeros(N); fullC=zeros(N); fullM=zeros(N);
for j=1:nx
 ids=(j-1)*n+(1:n); B=randn(stream,n,n); blocks{j}=B*B'/n+.3*eye(n);
 G{j}=randn(stream,n,n)/sqrt(n)+eye(n); M{j}=eye(n)+.08*randn(stream,n,n)/sqrt(n);
 fullG(ids,ids)=G{j}; fullC(ids,ids)=blocks{j}; fullM(ids,ids)=M{j};
end
fullC=fullC+U*U'; C=p1_covariance(blocks,U); z=randn(stream,N,2);
precisionError=norm(p1_precision(C,z)-fullC\z,'fro')/norm(fullC\z,'fro');
assert(precisionError<1e-11);
H=[speye(n),-speye(n);sparse(n,N)]; lambda=[.03 .04]; mu=.1; scale=.7;
opts=p1_options(struct('cov_floor',1e-6),struct('pcg_tolerance',1e-11,'pcg_max_iterations',300));
s0=randn(stream,N,1); b=.1*randn(stream,N,1); D=reshape(fullG*s0+b,n,nx);
Y=reshape(D(:)-b,n,nx); ctx=p1_context(G,C,Y,H,lambda,mu,opts);
[s,rec]=p1_solve(ctx,scale,zeros(n,nx));
ridge=repelem(lambda(:),n); A=fullG'*(fullC\fullG)+scale*diag(ridge)+scale*mu*(H'*H);
exact=A\(fullG'*(fullC\Y(:))); solutionError=norm(s(:)-exact)/norm(exact);
assert(solutionError<1e-8 && rec.normal_residual<2e-7);
stageBlocks=cell(1,nx); stageG=stageBlocks; stageU=fullM*U; stageY=reshape(fullM*Y(:),n,nx);
for j=1:nx, stageBlocks{j}=M{j}*blocks{j}*M{j}'; stageG{j}=M{j}*G{j}; end
Cs=p1_covariance(stageBlocks,stageU); sc=p1_context(stageG,Cs,stageY,H,lambda,mu,opts);
[ss,srec]=p1_solve(sc,scale,s); splitError=norm(ss(:)-s(:))/max(norm(s(:)),realmin);
assert(splitError<1e-8 && srec.normal_residual<2e-7);
floors=[.02 .03]; physical=fullM*fullC*fullM'; fullFloor=diag(repelem(floors,n));
for j=1:nx, stageBlocks{j}=stageBlocks{j}+floors(j)*eye(n); end
Cf=p1_covariance(stageBlocks,stageU,floors); exactTarget=trace((physical+fullFloor)\physical);
targetError=abs(Cf.target_squared-exactTarget); assert(targetError<1e-10);
% Noise-window bootstrap is deterministic and does not alter MATLAB global RNG.
t=(0:200)'*.002; noise=randn(stream,numel(t),3); po=p1_options(struct('cov_floor',1e-6),struct);
q1=p1_noise_summary(noise,t,po,47); q2=p1_noise_summary(noise,t,po,47);
assert(isequal(q1.bootstrap_variances,q2.bootstrap_variances));
assert(isequal(before,rng),'P1 modified the global MATLAB RNG');
assert(all(q1.variance_upper>=q1.variance_reference));
report=struct('status','PASS','precision_relative_error',precisionError, ...
 'quadratic_solution_relative_error',solutionError,'normal_residual',rec.normal_residual, ...
 'exact_propagated_split_relative_error',splitError,'split_normal_residual',srec.normal_residual, ...
 'floor_adjusted_target_absolute_error',targetError,'rng_preserved',true, ...
 'scope','small-matrix mathematics, fixed mean sign, full-gather Woodbury covariance, exact Split propagation and deterministic noise bootstrap; not geological validation');
if nargin>0 && ~isempty(outdir)
 assert(~isfolder(outdir),'Use a new selftest directory'); mkdir(outdir); uq_writejson(fullfile(outdir,'p1_selftest.json'),report);
end
fprintf('P1 small-matrix selftest PASS: precision %.3e, solution %.3e, Split %.3e\n',precisionError,solutionError,splitError);
end
