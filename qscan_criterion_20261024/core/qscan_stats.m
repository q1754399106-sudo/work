function out=qscan_stats(D,sigma,used,offset,grid)
% For every common factor c in grid: lasso statistics (same lambda rule as
% FullCalRefit/QScan) and the FullCalRefit output with Q*c, so that any
% selection criterion can be evaluated afterwards without re-solving.
[n,nx]=size(D);sigma=reshape(sigma,1,[]);live=any(D~=0,1);lambda=.5*sqrt(2*log(n))*sigma;nc=numel(grid);
J=zeros(nc,nx);RSSl=J;L1=J;RSSr=J;K=J;G=zeros(n,nx,nc,'single');flags=0;normal=0;
for q=1:nc
 c=used;c.Q=used.Q*grid(q);d=prune_dictionary(c,offset);rhs=d.B'*D(:,live);
 X=refit_lasso(d.Gram,rhs,lambda(live),d.L);res=D(:,live)-d.B*X;
 RSSl(q,live)=sum(res.^2,1);L1(q,live)=sum(abs(X),1);J(q,live)=.5*RSSl(q,live)+lambda(live).*L1(q,live);
 S=zeros(n,nx);
 for first=1:11:nx
  ids=first:min(first+10,nx);o=reflectivity_refit_prepared(D(:,ids),repmat({d},1,numel(ids)),sigma(ids));
  S(:,ids)=o.gather;K(q,ids)=cellfun(@numel,o.support);flags=flags+o.optimization_flags;normal=max(normal,o.normal_residual);
 end
 RSSr(q,:)=sum((D-d.G*S).^2,1);G(:,:,q)=single(S);
end
s2=sigma.^2;s2(~live)=Inf;
out=struct('grid',grid,'J_lasso',J./s2,'RSS_lasso',RSSl./s2,'L1_over_sigma',L1./sigma,'RSS_refit',RSSr./s2,'support',K, ...
 'live',live,'n',n,'gathers',G,'optimization_flags',flags,'max_normal_residual',normal);
end
