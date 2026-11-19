function fld_post()
% Post hoc check (PLAN_STAGE4_POSTHOC_ZH.md): pass 2 of fld_run.m with the dispersion reference
% frequency at 25 Hz (variant A) and amplitude only, k = 0 (variant B; f_ref = f). Same noise estimate,
% t* fields, methods and m. Outputs cache/fld/images_fref25.mat and images_amp.mat.
P=fld_setup();maxNumCompThreads(P.threads);ns=P.nshot;N=P.nz*P.nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([P.nz P.nx],Zi(:),Xi(:));nin=numel(lin);
irec=sub2ind([P.nz P.nx],P.rec_iz,P.rec_ix);isrc=sub2ind([P.nz P.nx],P.src_iz,P.src_ix);col=sub2ind([N ns],irec,P.shot_of);
S=load(fullfile(P.cache,'inv.mat'));INV=S.iv;TS=fld_tstar(P,lin);I0=load(fullfile(P.cache,'images.mat'),'meta');fcut=I0.meta.fcut;
meth={'PH','IQ30','IQ70','LP','W0','UQW0','W1','UQW1'};ms=[.5 .6 .7 .85 1 1.2 1.4 1.7 2];var={'fref25','amp'};
xmax=1.05*double(max(TS(:)))/min(ms);xg=single(linspace(0,xmax,1201));dxg=double(xg(2));s30=gain_s(30);s70=gain_s(70);
ck=fullfile(P.cache,'state_post.mat');Img=struct('fref25',struct(),'amp',struct());j0=0;
if isfile(ck),S=load(ck);Img=S.Img;j0=S.j;end
for j=j0+1:numel(P.f)
 t0=tic;f=P.f(j);[L,U,Pp,Qq,R]=lu(helm_op(P.vmig,f,P));
 b=zeros(N,ns);b(sub2ind([N ns],isrc,(1:ns)'))=-P.W(j)/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
 b=zeros(N,ns);b(col)=P.D(:,j);V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;clear L U Pp Qq R b
 iv=INV(:,j);Lf=single(f<=fcut-5)+single(f>fcut-5&f<fcut)*.5*(1+cos(pi*(f-fcut+5)/5));
 for v=1:2
  if v==1,fr=25;else,fr=f;end
  T0=uq_tab(f,xg,0,fr);T1=uq_tab(f,xg,P.sig,fr);
  acc=struct();for q=meth,acc.(q{1})=zeros(nin,numel(ms),'single');end
  for c=1:ceil(ns/31)
   sh=(c-1)*31+1:min(c*31,ns);Cs=Us(:,sh).*single(V(lin,sh));
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
  for q=meth,if isfield(Img.(var{v}),q{1}),Img.(var{v}).(q{1})=Img.(var{v}).(q{1})+acc.(q{1});else,Img.(var{v}).(q{1})=acc.(q{1});end,end
 end
 if mod(j,10)==0||j==numel(P.f),save(ck,'Img','j','-v7.3');end
 fprintf('fld post f %5.2f Hz (%d/%d): %.0f s\n',f,j,numel(P.f),toc(t0));
end
for v=1:2
 X=Img.(var{v});X.meta=struct('methods',{meth},'m',ms,'fcut',fcut,'variant',var{v});
 R0=load(fullfile(P.cache,'images.mat'),'R0');X.R0=R0.R0;save(fullfile(P.cache,['images_' var{v} '.mat']),'-struct','X','-v7.3');
end
end
