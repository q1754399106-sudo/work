% Verbatim from qcriteria_20261024/run_criteria.m (function criteria).
function v=criteria(D,A,d,sig)
% [C0 lasso objective; C1 Gaussian evidence; C2 GCV; C3 5-fold block hold-out]
[n,m]=size(D);lam=.5*sqrt(2*log(n))*sig;
X=refit_lasso(d.Gram,d.B'*D,lam,d.L);res=D-d.B*X;C0=sum((.5*sum(res.^2,1)+lam.*sum(abs(X),1))./sig.^2);
[U,S,~]=svd(A);s2=diag(S).^2;P=(U'*D).^2;s2grid=logspace(-8,2,81);C1=0;mu=zeros(1,m);
for j=1:m
 den=sig(j)^2+s2grid.*s2;ev=sum(log(den),1)+sum(P(:,j)./den,1);[e,i]=min(ev);C1=C1+.5*e;mu(j)=sig(j)^2/s2grid(i);
end
mugrid=logspace(-10,2,61)*max(s2);C2=0;
for j=1:m
 phi=s2./(s2+mugrid);r=sum(((1-phi).^2).*P(:,j),1);g=n*r./(n-sum(phi,1)).^2;C2=C2+min(g);
end
folds=round(linspace(0,n,6));C3=0;
for f=1:5
 te=folds(f)+1:folds(f+1);tr=setdiff(1:n,te);[Ut,St,Vt]=svd(A(tr,:),'econ');st=diag(St);B=Ut'*D(tr,:);
 for j=1:m,x=Vt*((st./(st.^2+mu(j))).*B(:,j));C3=C3+sum((A(te,:)*x-D(te,j)).^2)/sig(j)^2;end
end
v=[C0;C1;C2;C3];
end
