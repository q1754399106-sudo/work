function B=p1_layered_calibration_templates(public,offset,layer_q,fp)
% Nonstationary layered-Q basis. Travel coefficients depend only on public vp/h.
persistent cached_public cached_offsets cached_tau cached_E
assert(numel(layer_q)==4 && all(isfinite(layer_q)) && all(layer_q>0));
if isempty(cached_public)||~isequal(cached_public,public)
 cached_public=public;cached_offsets=[];cached_tau={};cached_E=fft(eye(public.nt),public.nfft,1);
end
index=find(cached_offsets==offset,1);
if isempty(index)
 c=public;c.Q=ones(1,4);[~,~,~,rayp]=ray_absorption(public.t,offset,c);
 bounds=[0,cumsum(2*public.h./public.vp(1:3)),Inf];
 tau=max(0,min(public.t,bounds(2:end))-bounds(1:end-1));
 tau(~isfinite(tau))=0;
 tau=tau./sqrt(max(1-(rayp.*reshape(public.vp,1,[])).^2,realmin));
 zero=public.t<=0;tau(zero,:)=0;tau(zero,1)=offset/public.vp(1);
 cached_offsets(end+1)=offset;cached_tau{end+1}=tau;index=numel(cached_offsets);
end
times=cumsum(2*public.h./public.vp(1:3));W=zeros(public.nt,3);
for k=1:3,z=pi*fp*(public.t-times(k));W(:,k)=(1-2*z.^2).*exp(-z.^2);end
active=any(abs(W)>1e-10,2);N=public.nfft;f=(0:N/2)'/(N*public.dt);
z=-pi*f+2i*f.*log(max(f,realmin)/public.fref);z(1)=0;
attenuation=cached_tau{index}(active,:)*(1./layer_q(:));
H=exp(z*attenuation');H(end,:)=real(H(end,:));H=[H;conj(H(end-1:-1:2,:))];
columns=real(ifft(H.*cached_E(:,active),[],1));Y=columns(1:public.nt,:)*W(active,:);
phase=[0;-1i*ones(N/2-1,1);0;1i*ones(N/2-1,1)];quadrature=real(ifft(fft(Y,N,1).*phase,[],1));
B=zeros(public.nt,6);B(:,1:2:end)=Y;B(:,2:2:end)=quadrature(1:public.nt,:);
end
