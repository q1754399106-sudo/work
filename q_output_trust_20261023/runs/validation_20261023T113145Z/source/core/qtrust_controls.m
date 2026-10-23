function c=qtrust_controls()
n=128;rng(20261024);v=randn(n,1);v=v/norm(v);B=eye(n);B(:,60)=v;B(:,61)=v;
d=struct('B',B,'W',B,'scale',ones(n,1),'frequency',30,'dt',.002);r=zeros(n,1);r([60,61])=[.6;.4];y=B*r;db=repmat({zeros(n)},1,4);Cq=.1^2*eye(4);
a=qtrust_refit(y,d,1e-4,r,[60,61],db,Cq,1,true);
assert(a.rank_after==1&&norm(a.gather-y)/norm(y)<1e-6&&a.max_protected_change_sd<=1+1e-7);
b=qtrust_refit(y,d,1e-4,r,[60,61],db,zeros(4),1,true);e=qtrust_refit(y,d,1e-4,r,[60,61],db,zeros(4),1,false);assert(norm(b.gather-e.gather)<1e-12);
z=qtrust_refit(zeros(n,1),d,0,zeros(n,1),[],db,Cq,1,true);assert(norm(z.gather)==0);
c=struct('duplicate_rank_after',a.rank_after,'duplicate_error',norm(a.gather-y)/norm(y),'zero_covariance_difference',norm(b.gather-e.gather),'zero_input',norm(z.gather)==0,'protected_change_sd',a.max_protected_change_sd,'passed',true,'scope','Oracle-support numerical controls, not performance evidence');
end
