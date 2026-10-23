function out=qtrust_refit(y,d,sigma,initialR,support,dB,Cq,kappa,useQ)
% Cumulative output-space trust constraints; no target/truth inputs.
n=numel(y);support=reshape(support,1,[]);m=numel(support);
if norm(y)==0||m==0
 out=struct('gather',zeros(n,1),'rank_before',m,'rank_after',0,'normal_residual',0,'objective_path',0,'protected_modes',0,'output_modes',0,'max_protected_change_sd',0,'guard_rejections',0,'reference_change_norm',0,'removed_modes',0,'kappa',kappa,'use_Q',useQ);return
end
sigma=max(sigma,1e-12*norm(y)/sqrt(n));C=d.B(:,support);H=d.W(:,support)./d.scale(support)';r0=initialR(:).*d.scale;
J=zeros(n,4);for k=1:4,J(:,k)=dB{k}*r0;end
if useQ
 [V,S]=eig((Cq+Cq')/2);E=J*V*diag(sqrt(max(diag(S),0)))/sigma;[U,S,~]=svd(E,'econ');s=diag(S);w=1./sqrt(1+s.^2)-1;
 Cw=C+U*(w.*(U'*C));yw=y+U*(w.*(U'*y));
else,Cw=C;yw=y;
end
gap=max(1,round(.5/(d.frequency*d.dt)));groups=cumsum([1,diff(support)>gap]);ng=max(groups);T=zeros(m);ranges=cell(1,ng);ranks=zeros(1,ng);next=1;
for g=1:ng
 idx=find(groups==g);[~,~,V]=svd(Cw(:,idx),'econ');k=numel(idx);dest=next:next+k-1;T(idx,dest)=V;ranges{g}=dest;ranks(g)=k;next=next+k;
end
Z=Cw*T;G=Z'*Z;G=(G+G')/2;b=Z'*yw;yty=sum(yw.^2);Ht=H*T;
[Uh,Sh,~]=svd(H,'econ');sh=diag(Sh);Uh=Uh(:,sh>max(sh)*1e-8);P=Uh'*Ht;
keep=1:m;[z,loss,nr,reg]=fit(keep);zref=z;sref=Ht*zref;
A=G+reg*eye(m);L=A\P';variance=max(0,sigma^2*sum(L.*(G*L),1)');sd=max(sqrt(variance),1e-8*sigma);
pref=P*zref;protected=abs(pref)>2*sd;penalty=sigma^2*log(n);path=loss+penalty*m;rejections=0;
while any(ranks>0)
 best=path(end);chosen=0;bestidx=[];bestz=[];bestnr=nr;
 for g=1:ng
  if ranks(g)==0,continue;end
  drop=ranges{g}(ranks(g));idx=keep(keep~=drop);[az,rss,normal]=fit(idx);objective=rss+penalty*numel(idx);
  if objective>=best-1e-12*max(yty,realmin),continue;end
  delta=P(:,idx)*az-pref;
  if any(abs(delta(protected))./sd(protected)>kappa+1e-9),rejections=rejections+1;continue;end
  best=objective;chosen=g;bestidx=idx;bestz=az;bestnr=normal;
 end
 if chosen==0,break;end
 ranks(chosen)=ranks(chosen)-1;keep=bestidx;z=bestz;nr=bestnr;path(end+1)=best;
end
coeff=zeros(m,1);coeff(keep)=z;s=Ht*coeff;delta=P*coeff-pref;
if any(protected),maxchange=max(abs(delta(protected))./sd(protected));else,maxchange=0;end
assert(maxchange<=kappa+1e-7&&all(isfinite(s))&&nr<1e-8);
assert(all(diff(path)<=1e-10*max(path(1),realmin)));
out=struct('gather',s,'rank_before',m,'rank_after',numel(keep),'normal_residual',nr,'objective_path',path, ...
 'protected_modes',sum(protected),'output_modes',numel(sd),'max_protected_change_sd',maxchange,'guard_rejections',rejections, ...
 'reference_change_norm',norm(s-sref)/max(norm(sref),realmin),'removed_modes',m-numel(keep),'kappa',kappa,'use_Q',useQ);
 function [a,rss,normal,ridge]=fit(idx)
  if isempty(idx),a=[];rss=yty;normal=0;ridge=0;return;end
  Q=G(idx,idx);rhs=b(idx);[V,S]=eig(Q);e=max(diag(S),0);top=max(e);bottom=min(e);ridge=max(1e-8*top,(top-1e4*bottom)/(1e4-1));
  a=V*((V'*rhs)./max(e+ridge,realmin));rss=max(0,yty-2*a'*rhs+a'*Q*a);normal=norm(Q*a+ridge*a-rhs)/max(norm(rhs),realmin);
 end
end
