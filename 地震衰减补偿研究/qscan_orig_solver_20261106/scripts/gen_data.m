function d = gen_data(L, Qtrue, seed)
% Data from the original forward operator: d = A(Qtrue, f0) r + noise; independent sparse random reflectivity per trace.
rng(seed,'twister');
r=zeros(L.n1,L.ntr);
for k=1:L.ntr
    idx=find(rand(L.n1,1)<0.03); idx=idx(idx>10 & idx<L.n1-10);
    r(idx,k)=randn(numel(idx),1);
end
[~,~,A]=qcdn_core(zeros(L.n1,1),L.dt,L.f0,1,1,Qtrue);
x=A*r;
d=x+L.noise*sqrt(mean(x(:).^2))*randn(size(x));
end
