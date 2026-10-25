function out=check_output_covariance_mc(dest)
% Independent dense-covariance calculation checks the API's conditional linear
% noise propagation. Known support/frozen pilot only; not coverage validation.
rng(20261025);n=64;m=12;sigma=.01;support=1:m;
[Q,~]=qr(randn(n,m),0);C=Q*diag(linspace(.6,1.4,m));[Uh,~]=qr(randn(n,m),0);H=Uh*diag(linspace(.7,1.5,m));
B=eye(n);B(:,support)=C;W=eye(n);W(:,support)=H;d=struct('B',B,'W',W,'scale',ones(n,1),'frequency',30,'dt',.002);
r=zeros(n,1);r(support)=linspace(1,2,m);y=C*r(support);J=.05*randn(n,4);Cq=.1^2*eye(4);db=cell(1,4);
for k=1:4,z=zeros(n);z(:,1)=J(:,k)/r(1);db{k}=z;end
Sigma0=sigma^2*eye(n)+J*Cq*J';[Vraw,Eraw]=eig(Sigma0/sigma^2);
raw_asymmetry=norm(Sigma0-Sigma0','fro')/norm(Sigma0,'fro');raw_orthogonality=norm(Vraw'*Vraw-eye(n),'fro');raw_imaginary=max(abs(imag(diag(Eraw))));
% Roundoff asymmetry in a matrix with repeated eigenvalues can make eig use
% its general nonsymmetric path. Enforce the exact symmetric covariance.
Sigma=(Sigma0+Sigma0')/2;[V,E]=eig(Sigma/sigma^2);F=V*diag(1./sqrt(diag(E)))*V';Cw=F*C;G=Cw'*Cw;G=(G+G')/2;
ev=eig(G);rho=max(1e-8*max(ev),(max(ev)-1e4*min(ev))/(1e4-1));L=H*((G+rho*eye(m))\(Cw'*F));theory=diag(Uh'*(L*Sigma*L')*Uh);
N=512;values=zeros(m,N);max_rank_change=0;
for j=1:N
 noise=sigma*randn(n,1)+J*chol(Cq,'lower')*randn(4,1);a=qtrust_refit(y+noise,d,sigma,r,support,db,Cq,0,true);
 assert(a.protected_modes==m&&a.rank_after==m,'Zero-budget case must remain at the linear reference');values(:,j)=Uh'*a.gather;max_rank_change=max(max_rank_change,m-a.rank_after);
end
empirical=var(values,0,2);relative=abs(empirical-theory)./theory;
out=struct('passed',median(relative)<.15&&max(relative)<.35,'draws',N,'median_relative_variance_error',median(relative),'maximum_relative_variance_error',max(relative), ...
 'analytic_variances',theory,'empirical_variances',empirical,'rank_change',max_rank_change,'raw_covariance_asymmetry',raw_asymmetry,'raw_eigenvector_orthogonality_error',raw_orthogonality,'raw_eigenvalue_imaginary_max',raw_imaginary, ...
 'scope','Conditional fixed-support/fixed-pilot Gaussian check using independently formed dense covariance; not seismic performance or uncertainty coverage after selection');
uq_writejson(dest,out);assert(out.passed);
end
