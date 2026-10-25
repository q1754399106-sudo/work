function d=wavelet_dictionary(cm,offset,phase_deg,ghost)
% refit_dictionary/prune_dictionary with a shift-invariant source wavelet:
% Ricker(cm.fpeak) rotated by phase_deg (same convention as the generator's
% phase_source) plus ghost terms ghost=[amp delay_s; ...]. phase 0, no ghost
% reproduces prune_dictionary exactly up to the wavelet construction.
t=cm.t(:);n=numel(t);dt=cm.dt;L=2*n;tc=(0:L-1)'*dt;c0=n*dt;
z=pi*cm.fpeak*(tc-c0);w=(1-2*z.^2).*exp(-z.^2);
if phase_deg~=0
 N=2^nextpow2(4*L);F=fft(w,N);F=F(1:N/2+1);rot=exp(-1i*deg2rad(phase_deg))*ones(size(F));rot(1)=1;rot(end)=cos(deg2rad(phase_deg));
 G=F.*rot;w=real(ifft([G;conj(G(end-1:-1:2))]));w=w(1:L);
end
if nargin>=4,base=w;for k=1:size(ghost,1),w=w+ghost(k,1)*interp1(tc,base,tc-ghost(k,2),'linear',0);end,end
W=zeros(n,n);
for i=1:n,W(:,i)=w(n+1-(i-1)+(0:n-1)-0);end   % column i: wavelet centred at t(i)
[a,b]=ray_absorption(t,offset,cm);Gm=attenuation_matrix(a+b,cm);B=Gm*W;
scale=sqrt(sum(B.^2,1));assert(all(scale>0));B=B./scale;Gram=B'*B;Gram=(Gram+Gram')/2;
d=struct('W',W,'G',Gm,'B',B,'Gram',Gram,'scale',scale(:),'L',norm(Gram,1), ...
 'key',sprintf('%.16g,',[reshape(cm.Q,1,[]),cm.fpeak,offset,phase_deg]),'frequency',cm.fpeak,'dt',dt);
end
