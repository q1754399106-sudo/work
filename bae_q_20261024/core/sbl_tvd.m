function out=sbl_tvd(D,d,sigma)
% SBL time-variant deconvolution baseline (Yuan et al. 2017, IEEE TGRS):
% hierarchical Gaussian prior on the Q-filtered dictionary, type-II ML via the
% fast marginal-likelihood sequential add/re-estimate/delete algorithm
% (Tipping & Faul 2003). Noise variance fixed to the shared estimate sigma^2.
[n,nx]=size(D);sigma=reshape(sigma,1,[]);Gram=d.Gram;S=zeros(n,nx);R=S;iters=zeros(1,nx);K=iters;
for j=1:nx
 if all(D(:,j)==0),continue;end
 beta=1/sigma(j)^2;Py=d.B'*D(:,j);diagG=diag(Gram);
 alpha=Inf(n,1);[~,i]=max(Py.^2./diagG);
 alpha(i)=diagG(i)/max(Py(i)^2/diagG(i)-1/beta,eps);
 for it=1:4000
  a=find(isfinite(alpha));Ga=Gram(:,a);
  Sigma=inv(diag(alpha(a))+beta*Gram(a,a));Sigma=(Sigma+Sigma')/2;mu=beta*Sigma*Py(a);
  T=Ga*Sigma;Sf=beta*diagG-beta^2*sum(T.*Ga,2);Qf=beta*Py-beta^2*(T*Py(a));
  s=Sf;q=Qf;in=isfinite(alpha);s(in)=alpha(in).*Sf(in)./(alpha(in)-Sf(in));q(in)=alpha(in).*Qf(in)./(alpha(in)-Sf(in));
  theta=q.^2-s;dL=-Inf(n,1);newalpha=alpha;
  re=in&theta>0;newalpha(re)=s(re).^2./theta(re);delta=1./newalpha(re)-1./alpha(re);
  dL(re)=Qf(re).^2./(Sf(re)+1./delta)-log(1+Sf(re).*delta);
  ad=~in&theta>0;newalpha(ad)=s(ad).^2./theta(ad);dL(ad)=(Qf(ad).^2-Sf(ad))./Sf(ad)+log(Sf(ad)./Qf(ad).^2);
  de=in&theta<=0;if nnz(in)==1,de(:)=false;end
  dL(de)=Qf(de).^2./(Sf(de)-alpha(de))-log(1-Sf(de)./alpha(de));
  [best,i]=max(dL);if ~(best>1e-6),break;end
  if de(i),alpha(i)=Inf;else,alpha(i)=newalpha(i);end
 end
 a=find(isfinite(alpha));Sigma=inv(diag(alpha(a))+beta*Gram(a,a));mu=beta*Sigma*Py(a);
 coeff=zeros(n,1);coeff(a)=mu;R(:,j)=coeff./d.scale;S(:,j)=d.W*R(:,j);iters(j)=it;K(j)=numel(a);
end
assert(all(isfinite(S(:))));
out=struct('gather',S,'reflectivity',R,'iterations',iters,'active',K,'converged',all(iters<4000));
end
