function out=reflectivity_refit_prepared(D,dicts,sigma)
% Same rule for every acquisition/model. No true event information.
[n,nx]=size(D);sigma=reshape(sigma,1,[]);lambda=.5*sqrt(2*log(n))*sigma;
X=zeros(n,nx);done=false(1,nx);lasso=cell(1,nx);S=zeros(n,nx);R=S;
supports=cell(1,nx);normal=zeros(1,nx);ridge=normal;before=normal;after=normal;
for j=1:nx
 if done(j),continue;end
 group=find(cellfun(@(d)strcmp(d.key,dicts{j}.key),dicts)&~done);d=dicts{j};rhs=d.B'*D(:,group);
 [x,di]=refit_lasso(d.Gram,rhs,lambda(group),d.L);X(:,group)=x;done(group)=true;
 for q=1:numel(group),jj=group(q);lasso{jj}=struct('iterations',di.iterations(q),'kkt_ratio',di.kkt_ratio(q),'coordinate_fallback',di.coordinate_fallback(q),'flag',di.kkt_ratio(q)>di.tolerance);end
end
for j=1:nx
 d=dicts{j};x=X(:,j);threshold=max(.001*max(abs(x)),.05*lambda(j));support=find(abs(x)>threshold);before(j)=numel(support);
 for pass=1:2
  if isempty(support),a=[];nr=0;reg=0;
  else
   C=d.B(:,support);[U,Z,V]=svd(C,'econ');z=diag(Z);top=max(z.^2);bottom=min(z.^2);
   reg=max(1e-8*top,(top-1e4*bottom)/(1e4-1));
   a=V*((z./(z.^2+reg)).*(U'*D(:,j)));
   nr=norm(C'*(C*a-D(:,j))+reg*a)/max(norm(C'*D(:,j)),realmin);
  end
  if pass==2,break;end
  residual=D(:,j);if ~isempty(support),residual=residual-d.B(:,support)*a;end
  score=abs(d.B'*residual)/max(sigma(j),realmin);
  peak=score>sqrt(2*log(n))&score>=[-Inf;score(1:end-1)]&score>=[score(2:end);-Inf];
  extra=setdiff(find(peak),support);if isempty(extra),break;end
  support=union(support,extra);
 end
 coeff=zeros(n,1);coeff(support)=a;R(:,j)=coeff./d.scale;S(:,j)=d.W*R(:,j);
 supports{j}=support;normal(j)=nr;ridge(j)=reg;after(j)=numel(support);
end
assert(all(isfinite(S(:)))&&max(normal)<1e-8);
flags=sum(cellfun(@(x)x.flag,lasso));
out=struct('gather',S,'reflectivity',R,'support',{supports},'support_before',before,'support_after',after, ...
 'lasso',{lasso},'optimization_flags',flags,'normal_residual',max(normal),'ridge',ridge,'lambda',lambda, ...
 'wavelet_assumption','zero-phase Ricker at calibrated frequency','development_only',true);
end
