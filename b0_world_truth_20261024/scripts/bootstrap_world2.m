function W=bootstrap_world2(X,Qgen,Qa,Qb,cm,offset,phi,sig,seeds)
% qscan_bias_diagnostic_20261024/scripts/bootstrap_world.m statement for statement, plus the
% mean output of the Qa estimator (W.mean_a), which serves as the bootstrap-mean reference.
c=cm;c.Q=Qgen;d=prune_dictionary(c,offset);Dh=d.G*X;K=numel(seeds);EX=sum(X(:).^2);
S1=zeros([size(X) 2]);S2=S1;E=zeros(K,2);
for k=1:K
 Dk=Dh+ar1_noise(phi,sig,size(X,1),seeds(k));sk=estimate_sigmas(Dk,cm.t);
 A=full_refit(Dk,sk,cm,Qa,offset);B=full_refit(Dk,sk,cm,Qb,offset);
 S1(:,:,1)=S1(:,:,1)+A;S1(:,:,2)=S1(:,:,2)+B;S2(:,:,1)=S2(:,:,1)+A.^2;S2(:,:,2)=S2(:,:,2)+B.^2;
 E(k,:)=[sum((A-X).^2,'all'),sum((B-X).^2,'all')];
end
m=S1/K;v=sum(max(S2/K-m.^2,0),[1 2])*K/(K-1);v=v(:)';b2=zeros(1,2);
for j=1:2,b2(j)=sum((m(:,:,j)-X).^2,'all')-v(j)/K;end
W=struct('mse',mean(E,1)/EX,'bias2',b2/EX,'var',v/EX,'rel_delta_mse',(mean(E(:,2))-mean(E(:,1)))/mean(E(:,1)),'mean_a',m(:,:,1));
end
