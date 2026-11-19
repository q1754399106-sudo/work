function fld_run()
% Stage 4 (PLAN_STAGE4_ZH.md): Volve line, frequency-domain two-way migration with the stage-3 operators.
% Pass 1: uncompensated and random-polarity single-frequency images -> band-summed 1/snr (as est_noise2).
% Pass 2: per-point operators for every method and Q multiplier m (as migrate2). Engine, eikonal, t*,
% operator tables, gain and noise estimate are the stage-3 functions; only geometry and data differ.
% Outputs cache/fld/inv.mat, cache/fld/images.mat (resumable state in cache/fld/state*.mat).
P=fld_setup();maxNumCompThreads(P.threads);out=P.cache;
ns=P.nshot;N=P.nz*P.nx;[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([P.nz P.nx],Zi(:),Xi(:));nin=numel(lin);
irec=sub2ind([P.nz P.nx],P.rec_iz,P.rec_ix);isrc=sub2ind([P.nz P.nx],P.src_iz,P.src_ix);
col=sub2ind([N ns],irec,P.shot_of);                                         % linear index into b (N x ns)
W=P.W;                                                                       % flat 4-25 Hz, cosine tails
rng(20261118,'twister');EPS=single(2*(rand(P.ntrace,1)>.5)-1);              % random polarity per trace
% ---------------- pass 1 ----------------
fi=fullfile(out,'inv.mat');
if ~isfile(fi)
 ck=fullfile(out,'state1.mat');I=zeros(nin,numel(P.f),'single');If=I;j0=0;
 if isfile(ck),S=load(ck);I=S.I;If=S.If;j0=S.j;end
 for j=j0+1:numel(P.f)
  t0=tic;f=P.f(j);[L,U,Pp,Qq,R]=lu(helm_op(P.vmig,f,P));
  b=zeros(N,ns);b(sub2ind([N ns],isrc,(1:ns)'))=-W(j)/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
  for q=1:2
   b=zeros(N,ns);d=P.D(:,j);if q==2,d=d.*EPS;end;b(col)=d;
   V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;img=real(sum(Us.*single(V(lin,:)),2));
   if q==1,I(:,j)=img;else,If(:,j)=img;end
  end
  if mod(j,10)==0||j==numel(P.f),save(ck,'I','If','j','-v7.3');end
  fprintf('fld pass 1 f %5.2f Hz (%d/%d): %.0f s\n',f,j,numel(P.f),toc(t0));
 end
 bw=round(400/P.dx);iv=zeros(nin,numel(P.f),'single');
 for j=1:numel(P.f)
  bnd=abs(P.f-P.f(j))<=5+1e-9;p=box(sum(I(:,bnd),2),P,bw);pf=box(sum(If(:,bnd),2),P,bw);
  iv(:,j)=min(pf./max(p-pf,realmin('single')),1e3);
 end
 R0=sum(I,2);save(fi,'iv','R0','-v7.3');clear I If
end
S=load(fi);INV=S.iv;
% ---------------- t* fields and two-way t* per shot ----------------
TS=fld_tstar(P,lin);
% LP cut-off: highest f of the run from 10 Hz with reservoir-window median snr >= 1
z=P.z(P.in{1});x=P.x(P.in{2});wr=(z>=2600&z<=3200)'&(x>=-2500&x<=1500);sm=1./median(INV(wr(:),:),1);
j0=find(P.f>=10,1);jj=j0;while jj<numel(P.f)&&sm(jj+1)>=1,jj=jj+1;end;if sm(j0)<1,jj=j0;end;fcut=P.f(jj);
fprintf('LP cut-off %.2f Hz\n',fcut);
% ---------------- pass 2 ----------------
meth={'PH','IQ30','IQ70','LP','W0','UQW0','W1','UQW1'};ms=[.5 .6 .7 .85 1 1.2 1.4 1.7 2];
xmax=1.05*double(max(TS(:)))/min(ms);xg=single(linspace(0,xmax,1201));dxg=double(xg(2));
s30=gain_s(30);s70=gain_s(70);ck=fullfile(out,'state2.mat');Img=struct();j0=0;
if isfile(ck),S=load(ck);Img=S.Img;j0=S.j;end
for j=j0+1:numel(P.f)
 t0=tic;f=P.f(j);[L,U,Pp,Qq,R]=lu(helm_op(P.vmig,f,P));
 b=zeros(N,ns);b(sub2ind([N ns],isrc,(1:ns)'))=-W(j)/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
 b=zeros(N,ns);b(col)=P.D(:,j);V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;clear L U Pp Qq R b
 T0=uq_tab(f,xg,0,P.fref);T1=uq_tab(f,xg,P.sig,P.fref);iv=INV(:,j);
 Lf=single(f<=fcut-5)+single(f>fcut-5&f<fcut)*.5*(1+cos(pi*(f-fcut+5)/5));
 acc=struct();acc.R0=zeros(nin,1,'single');for q=meth,acc.(q{1})=zeros(nin,numel(ms),'single');end
 for c=1:ceil(ns/31)                                                      % shot chunks (memory)
  sh=(c-1)*31+1:min(c*31,ns);Cs=Us(:,sh).*single(V(lin,sh));acc.R0=acc.R0+real(sum(Cs,2));
  for k=1:numel(ms)
   u=TS(:,sh)/ms(k)/dxg;i0=min(floor(u),numel(xg)-2);ww=u-single(i0);i0=i0+1;it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
   A0=it(single(T0.A));C0=it(single(T0.C));P0=it(single(T0.P));A1=it(single(T1.A));C1=it(single(T1.C));P1=it(single(T1.P));
   for q=meth
    switch q{1}
     case 'PH',H=P0;
     case 'IQ30',H=(A0+s30*P0.*C0)./(1+s30*C0);
     case 'IQ70',H=(A0+s70*P0.*C0)./(1+s70*C0);
     case 'LP',H=A0*Lf;
     case 'W0',H=A0./(1+s70*C0+iv);
     case 'UQW0',H=A1./(1+s70*C1+iv);
     case 'W1',H=(A0+s70*P0.*C0+P0.*iv)./(1+s70*C0+iv);
     case 'UQW1',H=(A1+s70*P1.*C1+P1.*iv)./(1+s70*C1+iv);
    end
    acc.(q{1})(:,k)=acc.(q{1})(:,k)+real(sum(H.*Cs,2));
   end
  end
 end
 for q=fieldnames(acc)',if isfield(Img,q{1}),Img.(q{1})=Img.(q{1})+acc.(q{1});else,Img.(q{1})=acc.(q{1});end,end
 if mod(j,10)==0||j==numel(P.f),save(ck,'Img','j','-v7.3');end
 fprintf('fld pass 2 f %5.2f Hz (%d/%d): %.0f s\n',f,j,numel(P.f),toc(t0));
end
Img.meta=struct('methods',{meth},'m',ms,'fcut',fcut);save(fullfile(out,'images.mat'),'-struct','Img','-v7.3');
end

function y=box(v,P,bw)
y=reshape(v,P.nzi,P.nxi).^2;y=movmean(movmean(y,bw,1),bw,2);y=y(:);
end
