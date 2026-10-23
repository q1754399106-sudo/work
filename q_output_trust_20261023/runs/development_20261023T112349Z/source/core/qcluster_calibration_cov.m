function out=qcluster_calibration_cov(D,offsets,p,est)
% Local profiled sensitivity, with explicit weak-scale regularization.
% This is not a calibrated posterior/CI. Frequency is a nuisance parameter.
if isfield(est,'layer_Q'),q=reshape(est.layer_Q,1,[]);else,q=repmat(est.effective_Q,1,4);end
fp=est.estimated_frequency;step=.02;
times=cumsum(2*p.h./p.vp(1:3));mask=false(p.nt,1);
for tm=times,mask=mask|abs(p.t-tm)<=.12;end
H=zeros(5);rss=0;ns=0;sigmas=zeros(1,numel(offsets));
for j=1:numel(offsets)
 ni=estimate_trace_snr(D(:,j),p.t);sigmas(j)=max(ni.sigma,realmin);
 B=p1_layered_calibration_templates(p,offsets(j),q,fp);X=B(mask,:)/sigmas(j);y=D(mask,j)/sigmas(j);
 a=(X'*X+1e-8*max(trace(X'*X)/6,realmin)*eye(6))\(X'*y);
 e=y-X*a;w=sqrt(min(1,1.5./max(abs(e),realmin)));Xw=w.*X;[U,S,~]=svd(Xw,'econ');s=diag(S);U=U(:,s>max(s)*1e-8);
 J=zeros(numel(y),5);
 for k=1:5
  qp=q;qm=q;fpplus=fp;fpminus=fp;
  if k<=4,qp(k)=q(k)*exp(step);qm(k)=q(k)*exp(-step);else,fpplus=fp*exp(step);fpminus=fp*exp(-step);end
  bp=p1_layered_calibration_templates(p,offsets(j),qp,fpplus);bm=p1_layered_calibration_templates(p,offsets(j),qm,fpminus);
  v=w.*((bp(mask,:)-bm(mask,:))*a/(2*step*sigmas(j)));J(:,k)=v-U*(U'*v);
 end
 H=H+J'*J;rss=rss+sum((w.*e).^2);ns=ns+numel(y)-numel(a);
end
scale=max(1,rss/max(ns-5,1));precision=H/scale+diag([ones(1,4)/.5^2,1/.2231435513^2]);
C=precision\eye(5);C=(C+C')/2;[V,S]=eig(C(1:4,1:4));raw=diag(S);limited=min(.5^2,max(.03^2,raw));Cq=V*diag(limited)*V';
out=struct('logQ_covariance',Cq,'raw_eigen_variances',raw,'used_eigen_variances',limited,'profile_hessian',H,'residual_inflation',scale,'calibration_sigma',sigmas,'frequency_marginalized_locally',true,'prior_logQ_std',.5,'scope','Regularized local profiled calibration sensitivity; not validated uncertainty coverage');
end
