function out=qcluster_refit(y,d,sigma,initialR,support,dB,Cq,method)
% Uses observed data, shared initial support and calibration uncertainty only.
n=numel(y);support=reshape(support,1,[]);m=numel(support);r0=initialR(:).*d.scale;
if norm(y)==0||m==0
 out=struct('gather',zeros(n,1),'rank_before',m,'rank_after',0,'clusters',0,'cluster_sizes',[],'cluster_ranks',[],'normal_residual',0,'objective_path',0,'model_error_to_noise',0,'coefficient_modes',0,'finite',true);return
end
sigma=max(sigma,1e-12*norm(y)/sqrt(n));C=d.B(:,support);J=zeros(n,4);
for k=1:4,J(:,k)=dB{k}*r0;end
useQ=any(strcmp(method,{'QAwareRidge','QAwareCluster'}));
if useQ
 [Vq,Sq]=eig((Cq+Cq')/2);E=J*Vq*diag(sqrt(max(diag(Sq),0)))/sigma;
 [Ue,Se,~]=svd(E,'econ');se=diag(Se);weights=1./sqrt(1+se.^2)-1;
 Cy=C+Ue*(weights.*(Ue'*C));yw=y+Ue*(weights.*(Ue'*y));
else,Cy=C;yw=y;se=zeros(4,1);
end
gap=max(1,round(.5/(d.frequency*d.dt)));groups=cumsum([1,diff(support)>gap]);ng=max(groups);
T=zeros(m);ranges=cell(1,ng);sizes=zeros(1,ng);svds=zeros(m,1);first=1;
for g=1:ng
 idx=find(groups==g);[~,S,V]=svd(Cy(:,idx),'econ');k=numel(idx);dest=first:first+k-1;
 T(idx,dest)=V;ranges{g}=dest;sizes(g)=k;svds(dest)=diag(S);first=first+k;
end
Z=Cy*T;G=Z'*Z;G=(G+G')/2;b=Z'*yw;yty=sum(yw.^2);ranks=sizes;
if strcmp(method,'GlobalTSVD')
 [U,S,V]=svd(C,'econ');ss=diag(S);keep=ss>=.01*max(ss);a=V(:,keep)*((U(:,keep)'*y)./ss(keep));
 residual=C*a-y;nr=norm(V(:,keep)'*(C'*residual))/max(norm(C'*y),realmin);path=sum(residual.^2);
 ranks=NaN(size(sizes));total=sum(keep);
else
 keep=1:m;[z,loss,nr]=fit(keep);penalty=sigma^2*log(n);path=loss+penalty*numel(keep);
 if any(strcmp(method,{'ClusterModesNoQ','QAwareCluster'}))
  while any(ranks>0)
   best=path(end);chosen=0;bestkeep=[];bestz=[];bestnr=nr;
   for g=1:ng
    if ranks(g)==0,continue;end
    drop=ranges{g}(ranks(g));trial=keep(keep~=drop);[az,rss,normres]=fit(trial);objective=rss+penalty*numel(trial);
    if objective<best-1e-12*max(sum(yw.^2),realmin),best=objective;chosen=g;bestkeep=trial;bestz=az;bestnr=normres;end
   end
   if chosen==0,break;end
   ranks(chosen)=ranks(chosen)-1;keep=bestkeep;z=bestz;nr=bestnr;path(end+1)=best;
  end
 end
 a=zeros(m,1);a(keep)=z;a=T*a;total=numel(keep);
end
coeff=zeros(n,1);coeff(support)=a;gather=d.W*(coeff./d.scale);
assert(all(isfinite(gather))&&nr<1e-8,'Nonfinite output or failed normal-equation check');
assert(all(diff(path)<=1e-10*max(path(1),realmin)),'Rank-selection objective increased');
out=struct('gather',gather,'rank_before',m,'rank_after',total,'clusters',ng,'cluster_sizes',sizes,'cluster_ranks',ranks,'normal_residual',nr,'objective_path',path,'model_error_to_noise',norm(se),'coefficient_modes',total,'finite',all(isfinite(gather)));
 function [a,loss,nr]=fit(idx)
  if isempty(idx),a=[];loss=yty;nr=0;return;end
  Q=G(idx,idx);rhs=b(idx);[V,S]=eig(Q);s=max(diag(S),0);top=max(s);bottom=min(s);
  reg=max(1e-8*top,(top-1e4*bottom)/(1e4-1));a=V*((V'*rhs)./max(s+reg,realmin));
  loss=max(0,yty-2*a'*rhs+a'*Q*a);nr=norm(Q*a+reg*a-rhs)/max(norm(rhs),realmin);
 end
end
