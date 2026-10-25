function metrics=score_structure(S,scene,c)
truth=scene.truth; metrics=struct;
metrics.RelativeError=norm(S-truth,'fro')/norm(truth,'fro');
ii=c.t>=.28 & c.t<=1.08;
metrics.SignalWindowError=norm(S(ii,:)-truth(ii,:),'fro')/norm(truth(ii,:),'fro');
nn=c.t>=.02 & c.t<=.25;
metrics.NoiseWindowRMS=sqrt(mean(S(nn,:).^2,'all'));
amp=zeros(3,size(S,2)); amptrue=amp; theta=amp; unstable=false(size(amp));
for j=1:size(S,2)
 W=scene.templates{j}; amp(:,j)=pinv(W)*S(:,j); amptrue(:,j)=scene.coeff{j}; theta(:,j)=scene.sin2{j};
 if cond(W'*W)>1e6, unstable(:,j)=true; end
end
valid=~unstable; metrics.EventAmplitudeError=norm(amp(valid)-amptrue(valid))/max(norm(amptrue(valid)),realmin);
metrics.UnidentifiableAmplitudes=sum(unstable,'all');
metrics.AVOInterceptError=NaN; metrics.AVOGradientError=NaN;
if endsWith(scene.name,'AVO')
 est=zeros(2,3); ref=est;
 for k=1:3
  B=[ones(size(S,2),1),theta(k,:)']; est(:,k)=B\amp(k,:)'; ref(:,k)=B\amptrue(k,:)';
 end
 metrics.AVOInterceptError=norm(est(1,:)-ref(1,:))/norm(ref(1,:));
 metrics.AVOGradientError=norm(est(2,:)-ref(2,:))/norm(ref(2,:));
end
metrics.WeakEventError=NaN; metrics.FaultNeighborhoodError=NaN; metrics.PinchoutTailError=NaN;
if strcmp(scene.name,'WeakSection')
 metrics.WeakEventError=norm(amp(2,:)-amptrue(2,:))/norm(amptrue(2,:));
elseif strcmp(scene.name,'FaultSection')
 k=5:8; metrics.FaultNeighborhoodError=norm(S(:,k)-truth(:,k),'fro')/norm(truth(:,k),'fro');
elseif strcmp(scene.name,'PinchoutSection')
 k=9:11; t=c.t>.60 & c.t<.71;
 metrics.PinchoutTailError=norm(S(t,k)-truth(t,k),'fro')/norm(truth,'fro');
end
end
