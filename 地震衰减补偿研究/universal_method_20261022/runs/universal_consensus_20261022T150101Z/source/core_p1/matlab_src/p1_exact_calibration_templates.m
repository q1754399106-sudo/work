function B=p1_exact_calibration_templates(public,offset,q,fp)
% Full nonstationary G(Q)*Ricker basis, matching the reconstruction operator.
% Only a bounded Ricker-tail truncation is used; no event-stationary shortcut.
persistent cached_public cached_offsets cached_travel cached_E
if isempty(cached_public)||~isequal(cached_public,public)
 cached_public=public;cached_offsets=[];cached_travel={};
 cached_E=fft(eye(public.nt),public.nfft,1);
end
index=find(cached_offsets==offset,1);
if isempty(index)
 c=public;c.Q=ones(1,4);[a,b]=ray_absorption(public.t,offset,c);
 cached_offsets(end+1)=offset;cached_travel{end+1}=a+b;index=numel(cached_offsets);
end
times=cumsum(2*public.h./public.vp(1:3));W=zeros(public.nt,3);
for k=1:3,z=pi*fp*(public.t-times(k));W(:,k)=(1-2*z.^2).*exp(-z.^2);end
active=any(abs(W)>1e-10,2);N=public.nfft;f=(0:N/2)'/(N*public.dt);
z=-pi*f+2i*f.*log(max(f,realmin)/public.fref);z(1)=0;
H=exp(z*cached_travel{index}(active)'/q);H(end,:)=real(H(end,:));H=[H;conj(H(end-1:-1:2,:))];
columns=real(ifft(H.*cached_E(:,active),[],1));Y=columns(1:public.nt,:)*W(active,:);
phase=[0;-1i*ones(N/2-1,1);0;1i*ones(N/2-1,1)];
quadrature=real(ifft(fft(Y,N,1).*phase,[],1));
B=zeros(public.nt,6);B(:,1:2:end)=Y;B(:,2:2:end)=quadrature(1:public.nt,:);
end
