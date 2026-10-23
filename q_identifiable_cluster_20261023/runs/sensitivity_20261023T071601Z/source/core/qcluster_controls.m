function checks=qcluster_controls()
% Mechanism checks use known supports on purpose, never counted as evaluation.
n=128;rng(20261023);v=randn(n,1);v=v/norm(v);B=eye(n);B(:,60)=v;B(:,61)=v;
d=struct('B',B,'W',B,'scale',ones(n,1),'frequency',30,'dt',.002);r=zeros(n,1);r([60,61])=[.6;.4];y=B*r;
db=repmat({zeros(n)},1,4);Cq=.1^2*eye(4);
a=qcluster_refit(y,d,1e-4,r,[60,61],db,Cq,'QAwareCluster');
assert(a.rank_after==1&&norm(a.gather-y)/norm(y)<1e-6,'Duplicate atoms must compress without deleting the waveform');
b=qcluster_refit(y,d,1e-4,r,[60,61],db,zeros(4),'ClusterModesNoQ');
c=qcluster_refit(y,d,1e-4,r,[60,61],db,zeros(4),'QAwareCluster');
assert(norm(b.gather-c.gather)<1e-12);
db{1}(:,60)=v;db{1}(:,61)=v;u=qcluster_refit(y,d,1e-4,r,[60,61],db,Cq,'QAwareCluster');
assert(u.finite&&u.rank_after<=u.rank_before&&u.model_error_to_noise>0);
checks=struct('duplicate_rank_before',2,'duplicate_rank_after',a.rank_after,'duplicate_waveform_error',norm(a.gather-y)/norm(y),'zero_Q_covariance_limit',norm(b.gather-c.gather),'positive_Q_uncertainty_finite',u.finite,'scope','Oracle-support numerical controls only, not performance evidence');
end
