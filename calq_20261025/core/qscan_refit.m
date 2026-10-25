function [S,info]=qscan_refit(D,sigma,used,offset)
% Data-driven common Q correction (literature route A): pick c minimising the
% gather-summed lasso objective J=(0.5||D-Bx||^2+lambda||x||_1)/sigma^2 under
% the same lambda rule as FullCalRefit, then run FullCalRefit with Q*c.
[n,nx]=size(D);sigma=reshape(sigma,1,[]);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;
grid=exp(linspace(log(.6),log(2.2),13));J=zeros(size(grid));
for q=1:numel(grid)
 c=used;c.Q=used.Q*grid(q);d=prune_dictionary(c,offset);rhs=d.B'*D(:,live);
 X=refit_lasso(d.Gram,rhs,lambda(live),d.L);res=D(:,live)-d.B*X;
 J(q)=sum((.5*sum(res.^2,1)+lambda(live).*sum(abs(X),1))./sigma(live).^2);
end
[~,best]=min(J);c=used;c.Q=used.Q*grid(best);d=prune_dictionary(c,offset);S=zeros(n,nx);flags=0;normal=0;
for first=1:11:nx
 ids=first:min(first+10,nx);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sigma(ids));
 S(:,ids)=o.gather;flags=flags+o.optimization_flags;normal=max(normal,o.normal_residual);
end
info=struct('grid',grid,'objective',J,'selected_factor',grid(best),'applied_Q',c.Q,'optimization_flags',flags,'max_normal_residual',normal);
end
