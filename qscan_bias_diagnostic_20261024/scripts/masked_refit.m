function [Dpred,S]=masked_refit(D,sigma,cm,Q,offset,obs)
% FullCalRefit fitted on the observed rows obs only: reflectivity_refit_prepared statement for
% statement (same lambda rule with n = nt, same threshold, residual-peak pass and ridge refit),
% with the dictionary rows restricted. Returns the predicted data on all rows and the output.
% With obs all true it reproduces full_refit exactly (checked in check_inputs).
c=cm;c.Q=Q;d=prune_dictionary(c,offset);B=d.B(obs,:);Gram=B'*B;Gram=(Gram+Gram')/2;L=norm(Gram,1);
[n,nx]=size(D);Do=D(obs,:);sigma=reshape(sigma,1,[]);lambda=.5*sqrt(2*log(n))*sigma;Dpred=zeros(n,nx);S=Dpred;
for first=1:11:nx
 ids=first:min(first+10,nx);X=refit_lasso(Gram,B'*Do(:,ids),lambda(ids),L);
 for q=1:numel(ids)
  j=ids(q);x=X(:,q);threshold=max(.001*max(abs(x)),.05*lambda(j));support=find(abs(x)>threshold);
  for pass=1:2
   if isempty(support),a=[];
   else
    C=B(:,support);[U,Z,V]=svd(C,'econ');z=diag(Z);top=max(z.^2);bottom=min(z.^2);
    reg=max(1e-8*top,(top-1e4*bottom)/(1e4-1));a=V*((z./(z.^2+reg)).*(U'*Do(:,j)));
   end
   if pass==2,break;end
   residual=Do(:,j);if ~isempty(support),residual=residual-B(:,support)*a;end
   score=abs(B'*residual)/max(sigma(j),realmin);
   peak=score>sqrt(2*log(n))&score>=[-Inf;score(1:end-1)]&score>=[score(2:end);-Inf];
   extra=setdiff(find(peak),support);if isempty(extra),break;end
   support=union(support,extra);
  end
  coeff=zeros(n,1);coeff(support)=a;Dpred(:,j)=d.B*coeff;S(:,j)=d.W*(coeff./d.scale);
 end
end
end
