function diag_snr2()
% Diagnostic (no compensation scores): band-summed random-polarity noise estimate.
% Frequencies 15-45 Hz every 1 Hz; datasets att (no noise) and m30 (cache data). Per frequency the
% uncompensated image I and the random-polarity image I_f (same signs as migrate.m) are stored;
% estimates use band sums over |f' - f| <= hb (2.5 and 5 Hz): snr = (<B^2> - <B_f^2>) / <B_f^2>
% in a 400 m box. Truth: band image of att (signal) against band image of m30 - att (noise).
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model('P1',P);ns=numel(P.ixs);nz=P.nz;nx=P.nx;N=nz*nx;
D=load(fullfile(P.cache,'data_P1.mat'),'att','m30');
js=find(P.f>=15&P.f<=45&mod(P.f,1)==0);nf=numel(js);
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([nz nx],Zi(:),Xi(:));nin=numel(lin);
irec=sub2ind([nz nx],P.izs*ones(size(P.ixr)),P.ixr);isrc=sub2ind([nz nx],P.izs*ones(size(P.ixs)),P.ixs);
rng(P.flip_seed,'twister');EPS=single(2*(rand(numel(P.ixr),ns)>.5)-1);
t=(1:P.fd.nt)'*P.fd.dt;a=(pi*P.fd.fpeak*(t-P.fd.t0)).^2;q=(1-2*a).*exp(-a);
I=zeros(nin,nf,4,'single');                                  % att, att flip, m30, m30 flip
for k=1:nf
 j=js(k);f=P.f(j);W=(1i*2*pi*f)*(q.'*exp(-1i*2*pi*t*f))*P.fd.dt;
 A=helm_op(M.vmig,f,P);[L,U,Pp,Qq,R]=lu(A);clear A
 b=zeros(N,ns);b(sub2ind([N ns],isrc,1:ns))=-W/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
 c=0;
 for ds={'att','m30'}
  Dd=squeeze(D.(ds{1})(:,j,:));
  for fl=[1 0]
   c=c+1;b=zeros(N,ns);if fl,b(irec,:)=Dd;else,b(irec,:)=Dd.*EPS;end
   V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;I(:,k,c)=real(sum(Us.*single(V(lin,:)),2));
  end
 end
 fprintf('f %.0f done\n',f);
end
z=P.z(P.in{1});x=P.x(P.in{2});wd=(z>=820&z<=1550)'&(x>=300&x<=2100);wl=(z>=700&z<=1100)'&(x>=600&x<=1200);bw=round(P.snr_box/P.dx);
f=P.f(js);
for hb=[2.5 5]
 fprintf('\nband half-width %.1f Hz: median over window of snr (true | est. on att | est. on m30)\n   f  |  deep: true   att    m30  |  lens: true   att    m30\n',hb);
 for k=1:nf
  bnd=abs(f-f(k))<=hb;S=sum(I(:,bnd,1),2);Nn=sum(I(:,bnd,3),2)-S;
  e1=est(sum(I(:,bnd,1),2),sum(I(:,bnd,2),2),P,bw);e3=est(sum(I(:,bnd,3),2),sum(I(:,bnd,4),2),P,bw);
  tr=box2(S.^2,P,bw)./box2(Nn.^2,P,bw);
  fprintf('%5.1f | %8.3g %8.3g %8.3g | %8.3g %8.3g %8.3g\n',f(k),median(tr(wd)),median(e1(wd)),median(e3(wd)),median(tr(wl)),median(e1(wl)),median(e3(wl)));
 end
end
save(fullfile(P.cache,'diag_snr2_P1.mat'),'I','js','-v7.3');
end

function s=est(b,bf,P,bw)
p=box2(b.^2,P,bw);pf=box2(bf.^2,P,bw);s=max(p-pf,0)./max(pf,realmin('single'));
end
function y=box2(v,P,bw)
y=reshape(v,P.nzi,P.nxi);y=movmean(movmean(y,bw,1),bw,2);y=y(:);
end
