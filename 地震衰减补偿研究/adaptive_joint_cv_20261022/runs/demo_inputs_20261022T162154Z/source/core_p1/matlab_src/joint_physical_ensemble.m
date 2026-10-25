function e=joint_physical_ensemble(cm,offsets)
% Equal-weight, layerwise log-Q scenarios. Calibrated scale is a proxy,
% not a posterior. Fixed identity data metric W=I, common signal coordinates.
u=min(.25,max(.05,cm.q_uncertainty));nl=numel(cm.Q);
patterns=[eye(nl);-eye(nl)];weights=ones(2*nl,1)/(2*nl);
e=cell(1,numel(offsets));cache=containers.Map('KeyType','double','ValueType','any');
for j=1:numel(offsets)
 x=offsets(j);if isKey(cache,x),e{j}=cache(x);continue;end
 Gs=cell(2*nl,1);Gbar=zeros(cm.nt);
 for k=1:2*nl
  ck=cm;ck.Q=cm.Q.*exp(sqrt(nl)*u*patterns(k,:));
  [a,b]=ray_absorption(cm.t,x,ck);Gs{k}=attenuation_matrix(a+b,ck);
  Gbar=Gbar+weights(k)*Gs{k};
 end
 R=zeros(cm.nt);for k=1:2*nl,dG=Gs{k}-Gbar;R=R+weights(k)*(dG'*dG);end
 R=(R+R')/2;
 probe=sin((1:cm.nt)'*.13);data=cos((1:cm.nt)'*.07);
 lhs=0;for k=1:2*nl,lhs=lhs+weights(k)*norm(Gs{k}*probe-data)^2;end
 rhs=norm(Gbar*probe-data)^2+probe'*R*probe;
 err=abs(lhs-rhs)/max(lhs,eps);assert(err<1e-11);
 e{j}=struct('Gbar',Gbar,'RQ',R,'identity_error',err,'logQ_scale',u,'weights',weights);
 cache(x)=e{j};
end
end
