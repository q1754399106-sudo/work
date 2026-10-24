function sel=qscan_select(D,sigma,cm,Q,offset)
% Selection step of bae_q_20261024/core/qscan_refit.m, statement for statement, with the two
% objective terms also kept. QScan output = full_refit at the selected factor (checked in stage 0).
[n,~]=size(D);sigma=reshape(sigma,1,[]);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;
used=cm;used.Q=Q;grid=exp(linspace(log(.6),log(2.2),13));J=zeros(size(grid));misfit=J;penalty=J;active=J;
for q=1:numel(grid)
 c=used;c.Q=used.Q*grid(q);d=prune_dictionary(c,offset);rhs=d.B'*D(:,live);
 X=refit_lasso(d.Gram,rhs,lambda(live),d.L);res=D(:,live)-d.B*X;
 J(q)=sum((.5*sum(res.^2,1)+lambda(live).*sum(abs(X),1))./sigma(live).^2);
 misfit(q)=sum(.5*sum(res.^2,1)./sigma(live).^2);penalty(q)=sum(lambda(live).*sum(abs(X),1)./sigma(live).^2);
 active(q)=nnz(X)/sum(live);
end
[~,best]=min(J);
sel=struct('grid',grid,'objective',J,'misfit',misfit,'penalty',penalty,'lasso_active',active, ...
 'best',best,'selected_factor',grid(best));
end
