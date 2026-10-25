function out=bae_q_refit(D,d,Gs,sigma,initial,use_mean)
% Bayesian approximation error (BAE) refit for uncertain Q.
% Model error eps_s=(G(Q_s)-G(Qhat))*s0, s0=FullCalRefit waveform (plug-in).
% Noise model sigma^2 I+Cov(eps); data shifted by mean(eps) when use_mean.
% Same lasso/threshold/peak-addition/ridge rules as reflectivity_refit_prepared,
% applied in the whitened, column-renormalised metric (unit noise).
[n,nx]=size(D);sigma=reshape(sigma,1,[]);lambda=.5*sqrt(2*log(n));M=numel(Gs);
S=zeros(n,nx);R=S;supports=cell(1,nx);normal=zeros(1,nx);flags=0;rank_used=normal;
for j=1:nx
 if all(D(:,j)==0),continue;end
 s0=initial.gather(:,j);E=zeros(n,M);
 for k=1:M,E(:,k)=Gs{k}*s0-d.G*s0;end
 mu=mean(E,2);if ~use_mean,mu=zeros(n,1);end
 [U,Z]=svd((E-mean(E,2))/sqrt(M-1),'econ');ev=diag(Z).^2;keep=ev>1e-12*max(max(ev),realmin);
 U=U(:,keep);shrink=1-1./sqrt(1+ev(keep)/sigma(j)^2);rank_used(j)=nnz(keep);
 white=@(v)(v-U*(shrink.*(U'*v)))/sigma(j);
 y=white(D(:,j)-mu);Bw=white(d.B);sw=sqrt(sum(Bw.^2,1));Bn=Bw./sw;
 Gram=Bn'*Bn;Gram=(Gram+Gram')/2;[x,di]=refit_lasso(Gram,Bn'*y,lambda,norm(Gram,1));flags=flags+di.flags;
 threshold=max(.001*max(abs(x)),.05*lambda);support=find(abs(x)>threshold);
 for pass=1:2
  if isempty(support),a=[];nr=0;
  else
   C=Bn(:,support);[Uc,Zc,Vc]=svd(C,'econ');z=diag(Zc);top=max(z.^2);bottom=min(z.^2);
   reg=max(1e-8*top,(top-1e4*bottom)/(1e4-1));a=Vc*((z./(z.^2+reg)).*(Uc'*y));
   nr=norm(C'*(C*a-y)+reg*a)/max(norm(C'*y),realmin);
  end
  if pass==2,break;end
  residual=y;if ~isempty(support),residual=residual-Bn(:,support)*a;end
  score=abs(Bn'*residual);
  peak=score>sqrt(2*log(n))&score>=[-Inf;score(1:end-1)]&score>=[score(2:end);-Inf];
  extra=setdiff(find(peak),support);if isempty(extra),break;end
  support=union(support,extra);
 end
 coeff=zeros(n,1);coeff(support)=a./sw(support)';R(:,j)=coeff./d.scale;S(:,j)=d.W*R(:,j);
 supports{j}=support;normal(j)=nr;
end
assert(all(isfinite(S(:)))&&max(normal)<1e-8);
out=struct('gather',S,'reflectivity',R,'support',{supports},'normal_residual',max(normal), ...
 'optimization_flags',flags,'error_rank',rank_used,'use_mean',use_mean, ...
 'assumption','Plug-in conditional BAE: model-error statistics from current waveform; Gaussian approximation; not a calibrated posterior.');
end
