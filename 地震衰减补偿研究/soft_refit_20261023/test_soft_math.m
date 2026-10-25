function test_soft_math()
% Verify propagated variance against direct multiplication and Monte Carlo.
rng(20261023);n=32;ids=[3;4;12;13;14];C=randn(n,numel(ids));C=C./sqrt(sum(C.^2,1));
W=randn(n,n)/sqrt(n);B=zeros(n,n);B(:,ids)=C;scale=ones(n,1);
dict=struct('B',B,'W',W,'scale',scale,'source_gram',W'*W,'frequency',25,'dt',.004);
y=C*[.4;-.1;.02;-.03;.015];[~,Z,~]=svd(C,'econ');z=diag(Z);
reg=max(1e-8*max(z.^2),(max(z.^2)-1e4*min(z.^2))/(1e4-1));
K=(C'*C+reg*eye(numel(ids)))\C';a=K*y;r=zeros(n,1);r(ids)=a;
initial=struct('support',{{ids}},'reflectivity',r,'gather',W*r,'optimization_flags',0);
sigma=.2;out=reflectivity_soft_prepared(y,{dict},sigma,initial);
group=cumsum([1;diff(ids)>round(.5/(25*.004))]);noise=sigma*randn(n,10000);
for k=1:max(group)
 ii=find(group==k);L=W(:,ids(ii))*K(ii,:);expected=sigma^2*sum(L(:).^2);
 reported=out.group_noise_energy{1}(k);
 assert(abs(reported-expected)<1e-10*max(expected,realmin));
 samples=L*noise;empirical=mean(sum(samples.^2,1));assert(abs(empirical/expected-1)<.08);
end
assert(all(abs(out.reflectivity(ids))<=abs(r(ids))));
assert(all(sign(out.reflectivity(ids))==sign(r(ids))));
zeroNoise=reflectivity_soft_prepared(y,{dict},0,initial);assert(norm(zeroNoise.gather-initial.gather)<1e-12);
empty=initial;empty.support={[]};empty.reflectivity=zeros(n,1);empty.gather=zeros(n,1);
dead=reflectivity_soft_prepared(zeros(n,1),{dict},realmin,empty);assert(all(dead.gather==0));
fprintf('Soft variance propagation, Monte Carlo, zero-noise and dead-trace checks passed.\n');
end
