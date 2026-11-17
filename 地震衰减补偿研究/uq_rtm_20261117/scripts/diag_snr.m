function diag_snr()
% Diagnostic of the image-domain SNR estimate (no compensation scores). Frequencies 20-40 Hz, every
% second one. Datasets: att (no noise) and att + noise at -30 dB (trace level, same definition as
% prep_data). Because migration is linear, the true noise image is I(noisy) - I(att).
% Compares, in W_deep and W_lens: true snr = <I_att^2>/<I_noise^2> per frequency and per 5 Hz band,
% against the even/odd estimate per frequency (stage 1) and per 5 Hz band sum.
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model('P1',P);ns=numel(P.ixs);nz=P.nz;nx=P.nx;N=nz*nx;
TD=load(fullfile(P.cache,'tdata_P1.mat'),'att');X=TD.att;rms_s=rmsdeep(X,P);
rng(777,'twister');n=randn(size(X,1),numel(X)/size(X,1));fq=(0:P.nt-1)'/(P.nt*P.dt);fq=min(fq,1/P.dt-fq);
n=real(ifft(fft(n).*(fq>=3&fq<=55)));n=reshape(single(n/sqrt(mean(n(:).^2))*rms_s/10^(-30/20)),size(X));
js=find(P.f>=20&P.f<=40);js=js(1:2:end);E=single(exp(-1i*2*pi*P.t*P.f(js)));
Da=dft(X,E);Dn=dft(X+n,E);
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([nz nx],Zi(:),Xi(:));
irec=sub2ind([nz nx],P.izs*ones(size(P.ixr)),P.ixr);isrc=sub2ind([nz nx],P.izs*ones(size(P.ixs)),P.ixs);W=wsrc(P,P.f(js));
z=P.z(P.in{1});x=P.x(P.in{2});wd=(z>=820&z<=1550)'&(x>=300&x<=2100);wl=(z>=700&z<=1100)'&(x>=600&x<=1200);
Ie=zeros(P.nzi*P.nxi,numel(js),2,'single');Io=Ie;
for q=1:numel(js)
 f=P.f(js(q));A=helm_op(M.vmig,f,P);[L,U,Pp,Qq,R]=lu(A);clear A
 b=zeros(N,ns);b(sub2ind([N ns],isrc,1:ns))=-W(q)/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
 for d=1:2
  if d==1,Dd=squeeze(Da(:,q,:));else,Dd=squeeze(Dn(:,q,:));end
  b=zeros(N,ns);b(irec,:)=Dd;V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;C=Us.*single(V(lin,:));
  Ie(:,q,d)=real(sum(C(:,1:2:end),2));Io(:,q,d)=real(sum(C(:,2:2:end),2));
 end
 fprintf('f %.1f done\n',f);
end
bw=round(P.snr_box/P.dx);
fprintf('\n   f   | true snr deep  lens | est.(per f) clean deep  noisy deep | est.(5 Hz band) clean deep  noisy deep | true band deep\n');
for q=1:numel(js)
 S=Ie(:,q,1)+Io(:,q,1);Nn=(Ie(:,q,2)+Io(:,q,2))-S;
 band=abs(P.f(js)-P.f(js(q)))<=2.5;
 Sb=sum(Ie(:,band,1)+Io(:,band,1),2);Nb=sum(Ie(:,band,2)+Io(:,band,2),2)-Sb;
 fprintf('%5.1f | %9.3g %9.3g | %9.3g %9.3g | %9.3g %9.3g | %9.3g\n',P.f(js(q)), ...
  mean(S(wd).^2)/mean(Nn(wd).^2),mean(S(wl).^2)/mean(Nn(wl).^2), ...
  med(est(Ie(:,q,1),Io(:,q,1),P,bw),wd),med(est(Ie(:,q,2),Io(:,q,2),P,bw),wd), ...
  med(est(sum(Ie(:,band,1),2),sum(Io(:,band,1),2),P,bw),wd),med(est(sum(Ie(:,band,2),2),sum(Io(:,band,2),2),P,bw),wd), ...
  mean(Sb(wd).^2)/mean(Nb(wd).^2));
end
end

function s=est(e,o,P,bw)
e=reshape(e,P.nzi,P.nxi);o=reshape(o,P.nzi,P.nxi);
s=4*max(movmean(movmean(e.*o,bw,1),bw,2),0)./max(movmean(movmean((e-o).^2,bw,1),bw,2),realmin('single'));
end
function y=med(s,w),y=median(s(w));end
function r=rmsdeep(X,P)
off=(P.ixr(:)-P.ixs(:)')*P.dx;tw=P.t>=1.2&P.t<=1.8;a=reshape(X(tw,:,:),nnz(tw),[]);r=sqrt(mean(a(:,abs(off(:))<=600).^2,'all'));
end
function Df=dft(X,E)
[nt,nr,ns]=size(X);Df=reshape(E.'*reshape(X,nt,nr*ns),[],nr,ns);Df=permute(Df,[2 1 3]);
end
function W=wsrc(P,f)
t=(1:P.fd.nt)'*P.fd.dt;a=(pi*P.fd.fpeak*(t-P.fd.t0)).^2;q=(1-2*a).*exp(-a);W=(1i*2*pi*f).*(q.'*exp(-1i*2*pi*t*f))*P.fd.dt;
end
