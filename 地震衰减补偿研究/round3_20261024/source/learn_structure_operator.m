function [H,info]=learn_structure_operator(D,c,o)
% Fixed local t-x prediction-error operator learned ONLY from observed D.
% AVO amplitudes are free: no sum-to-one condition on predictor coefficients.
% Correlation gating is an empirical edge safeguard, not a signal/noise oracle.
[n,nx]=size(D); f=(0:n-1)'/(n*c.dt); f=min(f,1/c.dt-f);
w=ones(n,1); low=f<o.pilot_pass_hz(1); high=f>o.pilot_pass_hz(2);
w(low)=0.5*(1-cos(pi*max(0,min(1,(f(low)-o.pilot_stop_hz(1))/ ...
 (o.pilot_pass_hz(1)-o.pilot_stop_hz(1))))));
w(high)=0.5*(1+cos(pi*max(0,min(1,(f(high)-o.pilot_pass_hz(2))/ ...
 (o.pilot_stop_hz(2)-o.pilot_pass_hz(2))))));
pilot=real(ifft(fft(D).*w)); lags=reshape(o.prediction_lags,1,[]);
centers=unique([1:o.prediction_step_samples:n,n]);
row=[]; col=[]; val=[]; r2all=[]; coeffall=[]; weights=zeros(n,nx);
half=floor(o.prediction_window_samples/2);
for j=1:nx
 neighbors=[j-1,j+1]; neighbors=neighbors(neighbors>=1 & neighbors<=nx);
 if isempty(neighbors), continue; end
 p=numel(neighbors)*numel(lags); B=zeros(n,p); ids=zeros(n,p); k=0;
 for jj=neighbors
  for lag=lags
   k=k+1; tt=min(n,max(1,(1:n)'+lag)); B(:,k)=pilot(tt,jj); ids(:,k)=tt+(jj-1)*n;
  end
 end
 C=zeros(numel(centers),p); gate=zeros(numel(centers),1);
 for k=1:numel(centers)
  ii=max(1,centers(k)-half):min(n,centers(k)+half); Z=B(ii,:); y=pilot(ii,j);
  ridge=o.prediction_ridge*max(sum(Z(:).^2)/p,realmin);
  coef=(Z'*Z+ridge*eye(p))\(Z'*y);
  coef=coef*min(1,o.prediction_coefficient_norm_cap/max(norm(coef),realmin));
  r2=1-sum((y-Z*coef).^2)/max(sum(y.^2),realmin);
  C(k,:)=coef'; gate(k)=sqrt(max(0,min(1,(r2-o.prediction_r2_gate)/(1-o.prediction_r2_gate))));
  r2all(end+1)=r2; coeffall(end+1)=norm(coef);
 end
 coef=interp1(centers,C,(1:n)','linear'); g=interp1(centers,gate,(1:n)','linear');
 g=g./sqrt(1+sum(coef.^2,2)); weights(:,j)=g;
 rr=(1:n)'+(j-1)*n; row=[row;rr]; col=[col;rr]; val=[val;g];
 for k=1:p
  row=[row;rr]; col=[col;ids(:,k)]; val=[val;-g.*coef(:,k)];
 end
end
H=sparse(row,col,val,n*nx,n*nx);
info=struct('mean_weight',mean(weights(:)),'active_fraction',mean(weights(:)>0), ...
 'r2_mean',mean(r2all),'max_coefficient_norm',max(coeffall),'weights',weights, ...
 'pilot',pilot,'nnz',nnz(H));
end
