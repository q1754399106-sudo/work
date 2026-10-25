function B=uq_templates(p,offset,q,fp)
% Event-local stationary propagation approximation, not G applied columnwise.
c=p; c.Q=ones(1,4); times=cumsum(2*p.h./p.vp(1:3));
[a,b]=ray_absorption(times,offset,c); travel=a+b;
N=p.nfft; f=(0:N/2)'/(N*p.dt); z=-pi*f+2i*f.*log(max(f,realmin)/p.fref); z(1)=0;
B=zeros(p.nt,6);
for k=1:3
 x=pi*fp*(p.t-times(k)); w=(1-2*x.^2).*exp(-x.^2);
 wf=fft(w,N); h=exp(z*travel(k)/q); h(end)=real(h(end)); h=[h;conj(h(end-1:-1:2))];
 u=wf.*h; v=ifft(u); B(:,2*k-1)=real(v(1:p.nt));
 phase=[0;-1i*ones(N/2-1,1);0;1i*ones(N/2-1,1)]; v=ifft(u.*phase); B(:,2*k)=real(v(1:p.nt));
end
end
