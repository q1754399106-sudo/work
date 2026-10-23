function [X,diagnostic]=refit_lasso(Gram,rhs,lambda,L)
% FISTA with a data-independent stopping tolerance; coordinate fallback.
[n,nx]=size(rhs);lambda=reshape(lambda,1,[]);X=zeros(n,nx);Y=X;theta=1;
ratio=Inf(1,nx);iterations=zeros(1,nx);tolerance=.02;
for k=1:3000
 previous=X;gradient=Gram*Y-rhs;Z=Y-gradient/L;X=sign(Z).*max(abs(Z)-lambda/L,0);
 next=(1+sqrt(1+4*theta^2))/2;Y=X+((theta-1)/next)*(X-previous);theta=next;
 if mod(k,25)==0
  ratio=kkt(Gram,X,rhs,lambda);finished=ratio<=tolerance;iterations(finished&iterations==0)=k;
  if all(finished),break;end
 end
end
fallback=false(1,nx);
for j=1:nx
 if ratio(j)<=tolerance,continue;end
 fallback(j)=true;x=X(:,j);lam=lambda(j);b=rhs(:,j);
 for outer=1:20
  gradient=Gram*x-b;active=find(abs(x)>0|abs(gradient)>lam*(1+tolerance/2));
  if isempty(active),break;end
  A=Gram(active,active);v=x(active);u=b(active)-A*v;
  for sweep=1:1200
   change=0;
   for q=1:numel(active)
    z=v(q)+u(q)/A(q,q);value=sign(z)*max(abs(z)-lam/A(q,q),0);delta=value-v(q);
    if delta~=0,v(q)=value;u=u-A(:,q)*delta;change=max(change,abs(delta));end
   end
   if change<1e-6*max(lam,max(abs(v))),break;end
  end
  x(:)=0;x(active)=v;
  if kkt(Gram,x,b,lam)<=tolerance,break;end
 end
 X(:,j)=x;
end
ratio=kkt(Gram,X,rhs,lambda);iterations(iterations==0)=k;
diagnostic=struct('iterations',iterations,'kkt_ratio',ratio,'coordinate_fallback',fallback,'flags',sum(ratio>tolerance),'tolerance',tolerance);
end
function ratio=kkt(A,X,b,lambda)
g=A*X-b;pen=ones(size(X,1),1)*reshape(lambda,1,[]);e=max(abs(g)-pen,0);active=X~=0;
e(active)=abs(g(active)+pen(active).*sign(X(active)));ratio=max(e,[],1)./reshape(lambda,1,[]);
end
