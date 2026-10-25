function [s,lambda,residual] = discrepancy_solve(G,y,grid,target)
% Select residual closest to supplied expected noise norm, without truth.
[U,S,V]=svd(G,'econ'); z=diag(S); uy=U'*y;
residual=zeros(size(grid));
for k=1:numel(grid)
    sk=V*((z./(z.^2+grid(k))).*uy);
    residual(k)=norm(G*sk-y);
end
[~,k]=min(abs(residual-target)); lambda=grid(k);
s=V*((z./(z.^2+lambda)).*uy);
end
