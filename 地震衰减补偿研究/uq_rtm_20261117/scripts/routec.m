function routec(name)
% Route C: data-domain UQ-W (gamma 0, sig P.sig, 70 dB floor) on every prestack trace, then the
% same acoustic migration as route A without any further operator. Two-way t*(t) of a trace:
% along the column below the source-receiver midpoint, source leg + receiver leg (flat-reflector
% specular point), from the same eikonal fields as route A. 1/snr(t, f): coherence of adjacent
% traces of the zero-offset section (0.4 s windows every 0.1 s, as UQ-W), one function for the line.
% Adds Df fields C_mXXX (m = Q multiplier x 100) to cache/data_<name>.mat.
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model(name,P);F=tstar_fields(M,P);
TD=load(fullfile(P.cache,['tdata_' name '.mat']),P.routeC.data);X=TD.(P.routeC.data);
ns=numel(P.ixs);nr=numel(P.ixr);N=1024;fo=(0:N/2)'/(N*P.dt);fo(1)=1e-6;n=P.nt;
% 1/snr on the zero-offset section
Z=zeros(n,ns);for k=1:ns,Z(:,k)=X(:,P.ixr==P.ixs(k),k);end
starts=P.t(1):P.snr_hop:P.t(end)-P.snr_win+1e-9;cen=starts+P.snr_win/2;inv=zeros(numel(starts),numel(fo));
for w=1:numel(starts)
 [f1,snr]=coherence_snr(Z,P.t,starts(w),starts(w)+P.snr_win,P.dt);inv(w,:)=interp1(f1,1./max(snr,P.snr_min),fo,'linear','extrap');
end
INV=interp1(cen(:),inv,min(max(P.t,cen(1)),cen(end)));      % n x nf
% operator tables over two-way t*
xmax=1.05*double(max(F.ts(:)))*2/min(P.routeC.m);xg=linspace(0,xmax,801);dxg=xg(2);
TA=zeros(numel(fo),numel(xg));TC=TA;
for j=1:numel(fo),T=uq_tab(fo(j),xg,P.sig,P.fd.fref);TA(j,:)=T.A;TC(j,:)=T.C;end
s=gain_s(70);E2=exp(1i*2*pi*(0:n-1)'*fo');cw=[1;2*ones(N/2-1,1);1]/N;
D=load(fullfile(P.cache,['data_' name '.mat']));Eh=single(exp(-1i*2*pi*P.t*P.f));
xs=P.x(P.ixs);rows=P.in{1};
for m=P.routeC.m
 tag=sprintf('C_m%03d',round(100*m));if isfield(D,tag),continue;end
 Y=zeros(n,nr,ns,'single');t0=tic;
 for k=1:ns
  for r=1:nr
   xr=P.x(P.ixr(r));jm=round(((xs(k)+xr)/2)/P.dx)+P.npml+1;
   q=min(max((xr-xs(1))/(xs(2)-xs(1))+1,1),ns);k0=min(floor(q),ns-1);wq=q-k0;
   Tc=double(F.T(rows,jm,k)+(1-wq)*F.T(rows,jm,k0)+wq*F.T(rows,jm,k0+1));
   tc=double(F.ts(rows,jm,k)+(1-wq)*F.ts(rows,jm,k0)+wq*F.ts(rows,jm,k0+1));
   Tc=cummax(Tc);[Tc,iu]=unique(Tc);xh=max(interp1(Tc,tc(iu),P.t,'linear','extrap'),0)/m;
   u=xh/dxg;i0=min(floor(u),numel(xg)-2);ww=min(u-i0,1);i0=i0+1;
   A=TA(:,i0).'.*(1-ww)+TA(:,i0+1).'.*ww;C=TC(:,i0).'.*(1-ww)+TC(:,i0+1).'.*ww;
   H=A./(1+s*C+INV);H(:,1)=real(H(:,1));H(:,end)=real(H(:,end));
   d=fft(double(X(:,r,k)),N);d=d(1:N/2+1);
   Y(:,r,k)=real((H.*E2)*(cw.*d));
  end
 end
 Df=reshape(Eh.'*reshape(Y,n,[]),[],nr,ns);D.(tag)=permute(Df,[2 1 3]);
 S=struct(tag,D.(tag));save(fullfile(P.cache,['data_' name '.mat']),'-struct','S','-append');
 fprintf('route C %s %s: %.0f s\n',name,tag,toc(t0));
end
end

function [f,snr]=coherence_snr(D,t,t0,t1,dt)
% as field_application_20261028/scripts/lib.py coherence_snr (nfft 512, Tukey 0.2, 5-point smoothing)
w=t>=t0-1e-9&t<=t1+1e-9;n=nnz(w);m=floor(.2*(n-1)/2);tp=ones(n,1);
if m>0,r=.5*(1-cos(pi*(0:m)'/m));tp(1:m+1)=r;tp(end-m:end)=flipud(r);end
Xf=fft(D(w,:).*tp,512);Xf=Xf(1:257,:);f=(0:256)'/(512*dt);
Pw=movmean(mean(abs(Xf).^2,2),5);C=movmean(mean(real(Xf(:,1:end-1).*conj(Xf(:,2:end))),2),5);
snr=max(C,0)./max(Pw-max(C,0),1e-30*max(Pw));
end
