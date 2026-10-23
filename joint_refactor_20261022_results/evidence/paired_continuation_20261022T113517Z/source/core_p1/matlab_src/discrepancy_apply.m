function [s,lambda,residual,endpoint] = discrepancy_apply(F,Y,grid,target)
% Batched equivalent of discrepancy_solve. Never uses clean signal for tuning.
% Orthogonal residual is retained for rectangular forward operators.
u=F.U'*Y; perpendicular=max(0,sum(Y.^2,1)-sum(u.^2,1));
nr=size(Y,2); rr=zeros(numel(grid),nr);
for k=1:numel(grid)
    rr(k,:)=sqrt(sum((grid(k)./(F.z.^2+grid(k)).*u).^2,1)+perpendicular);
end
[~,ix]=min(abs(rr-target),[],1); lambda=reshape(grid(ix),1,[]);
s=F.V*((F.z./(F.z.^2+lambda)).*u);
residual=rr(sub2ind(size(rr),ix,1:nr));
endpoint=(ix==1 | ix==numel(grid));
end
