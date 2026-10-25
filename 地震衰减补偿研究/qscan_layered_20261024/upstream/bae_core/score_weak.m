function [err,gain]=score_weak(A,reflectivity,events,weak_columns,fpeak,dt,traces)
% Ricker-template projection at registered weak events (copied from coherence_soft core).
if nargin<7,traces=1:size(A,2);end
est=[];truth=[];n=size(A,1);
for j=reshape(traces,1,[])
 for e=reshape(weak_columns,1,[])
  k=events(j,e);if ~isfinite(k)||reflectivity(k,j)==0,continue;end
  ii=max(1,k-10):min(n,k+10);z=pi*fpeak*(ii'-k)*dt;w=(1-2*z.^2).*exp(-z.^2);
  est(end+1,1)=(A(ii,j)'*w)/(w'*w);truth(end+1,1)=reflectivity(k,j); %#ok<AGROW>
 end
end
err=norm(est-truth)/norm(truth);gain=(est'*truth)/(truth'*truth);
end
