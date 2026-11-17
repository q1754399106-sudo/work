function test_engine()
% Engine checks on the homogeneous model T0 only (no scores):
%  E1 Helmholtz Green's function x W against the elastic FD direct wave (one complex scale fitted)
%  E2 sign and size of attenuation and dispersion: FD visco / elastic against exp[(-pi f + i k) t*]
%  E3 eikonal t* against distance/(v Q)
%  E4 uq_tab: sig = 0, gamma = 1, 1/snr = 0 equals exp(-ikx)(b+s)/(b^2+s); sig > 0 finite
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model('T0',P);
fd_shots('T0');S=load(fullfile(P.cache,'shots_T0','s000.mat'));
ixs=P.npml+3;cols=P.in{2};off=(cols-ixs)*P.dx;use=off>=40&off<=2200;
t=(1:P.nt)'*P.dt;E=exp(-1i*2*pi*t*P.f);                      % nt x nf
De=(double(S.ela).'*E);Da=(double(S.att).'*E);              % nrec x nf
W=wsrc(P);
% E1
jj=1:8:numel(P.f);Hd=zeros(numel(cols),numel(P.f));
for j=jj
 A=helm_op(M.v,P.f(j),P);b=zeros(P.nz*P.nx,1);b(sub2ind([P.nz P.nx],P.izs,ixs))=-W(j)/P.dx^2;
 u=A\b;Hd(:,j)=u(sub2ind([P.nz P.nx],P.izs*ones(size(cols)),cols));
end
a=Hd(use,jj);d=De(use,jj);al=(a(:)'*d(:))/(a(:)'*a(:));r=norm(al*a(:)-d(:))/norm(d(:));
fprintf('E1 Helmholtz vs FD elastic: relative misfit %.3f (scale %.3g%+.3gi)\n',r,real(al),imag(al));
for o=[200 1000 2000]
 i=find(off==o);x=al*Hd(i,jj);y=De(i,jj);
 fprintf('   offset %4d m: misfit %.3f\n     f   %s\n     dph %s\n     amp %s\n',o,norm(x-y)/norm(y), ...
  sprintf('%6.1f',P.f(jj)),sprintf('%6.2f',angle(x.*conj(y))),sprintf('%6.2f',abs(x)./abs(y)));
end
% E1b: the same at depth 800 m (away from the absorbing strip), receivers every 40 m
izd=P.npml+1+200;cd_=P.in{2}(1:10:end);q=(1-2*(pi*P.fd.fpeak*((1:P.fd.nt)'*P.fd.dt-P.fd.t0)).^2).*exp(-(pi*P.fd.fpeak*((1:P.fd.nt)'*P.fd.dt-P.fd.t0)).^2);
opt=struct('dx',P.dx,'dt',P.fd.dt,'nt',P.fd.nt,'fref',P.fd.fref,'npml',P.npml,'src',[P.izs ixs],'rec',[izd*ones(numel(cd_),1) cd_(:)]);
Me=M;Me.Q(:)=Inf;rd=fd2d_visco(Me,q,opt);rd=rd(10:10:end,:);Dd=rd.'*E;Hb=zeros(numel(cd_),numel(jj));
for n=1:numel(jj)
 j=jj(n);A=helm_op(M.v,P.f(j),P);b=zeros(P.nz*P.nx,1);b(sub2ind([P.nz P.nx],P.izs,ixs))=-W(j)/P.dx^2;
 u=A\b;Hb(:,n)=u(sub2ind([P.nz P.nx],izd*ones(size(cd_)),cd_));
end
x=al*Hb;y=Dd(:,jj);
fprintf('E1b depth 800 m, x 0-2360 m (same scale): relative misfit %.3f; per receiver max %.3f\n',norm(x(:)-y(:))/norm(y(:)),max(vecnorm(x-y,2,2)./vecnorm(y,2,2)));
% E2
R=Da(use,:)./De(use,:);ts=off(use)'/(1850*80);fb=P.f>=8&P.f<=50;
k=2*P.f.*log(P.f/P.fd.fref);model=exp((-pi*P.f+1i*k).*ts);
amp=-log(abs(R(:,fb)))./(pi*P.f(fb));ph=angle(R(:,fb).*conj(exp(1i*k(fb).*ts)));
fprintf('E2 fitted t*/model t*: median %.3f (5-95%%: %.3f-%.3f); residual phase |median| %.3f rad\n', ...
 median(amp./ts,'all'),prctile(amp./ts,5,'all'),prctile(amp./ts,95,'all'),abs(median(ph(:))));
fprintf('   |FD ratio - model| / |model| (8-50 Hz): median %.3f\n',median(abs(R(:,fb)-model(:,fb))./abs(model(:,fb)),'all'));
% E3
s=1./M.vmig;Tp=eikonal_fs(s.*(1+.5./M.Qm),P.dx,P.izs,ixs);Tm=eikonal_fs(s.*(1-.5./M.Qm),P.dx,P.izs,ixs);
T=eikonal_fs(s,P.dx,P.izs,ixs);tse=(Tp-Tm);[X,Z]=meshgrid(P.x,P.z);dist=hypot(X-P.x(ixs),Z-P.z(P.izs));
in=dist>100&dist<1500;
fprintf('E3 eikonal T/(dist/v) median %.4f; t*/(dist/(vQ)) median %.4f, max |rel err| %.4f\n', ...
 median(T(in)./(dist(in)/1850)),median(tse(in)./(dist(in)/(1850*80))),max(abs(tse(in)./(dist(in)/(1850*80))-1)));
% E4
xg=linspace(0,.1,401);f=40;t0=uq_tab(f,xg,0,P.fd.fref);s=gain_s(70);
H=(t0.A+s*t0.P.*t0.C)./(1+s*t0.C);b=exp(-pi*f*xg);Href=exp(-1i*2*f*log(f/P.fd.fref)*xg).*(b+s)./(b.^2+s);
t1=uq_tab(f,xg,P.sig,P.fd.fref);
fprintf('E4 sig=0 operator vs stabilised inverse Q: max rel %.1e; sig=%.3f finite %d; |H| max %.1f dB at 70 dB\n', ...
 max(abs(H-Href))/max(abs(Href)),P.sig,all(isfinite([t1.A t1.C t1.P])),20*log10(max(abs(H))));
end

function W=wsrc(P)
% spectrum of the wave-equation source of fd2d_visco: d/dt of the injected Ricker (fft convention)
t=(1:P.fd.nt)'*P.fd.dt;q=(1-2*(pi*P.fd.fpeak*(t-P.fd.t0)).^2).*exp(-(pi*P.fd.fpeak*(t-P.fd.t0)).^2);
W=(1i*2*pi*P.f).*(q.'*exp(-1i*2*pi*t*P.f))*P.fd.dt;
end
