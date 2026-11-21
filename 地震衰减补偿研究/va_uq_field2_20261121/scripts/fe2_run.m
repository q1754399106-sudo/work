function fe2_run(mode)
% Second field pass (PLAN_ZH.md section 5). Stage-4 engine (fld_setup, helm_op, fld_tstar, uq_tab, gain_s)
% and the step-1 VA-UQ table (va_tab), all unchanged; only the noise estimate (B = 300) and the methods
% are new. One MATLAB process, one LU per frequency, resumable.
%  pass 1, per frequency: LU of the stage-4 operator, source wavefield, back-propagation of
%    (a) the data: image products Cs = Us.*V (nin x nshot) written to cache/cs/fNNN.bin, image I(:,j);
%    (b) the block-flipped data (est_block.m rule, B = 300 consecutive receivers of a shot share a sign,
%        rng(20261118 + 300)): flipped image If(:,j);
%    then 1/snr exactly as fld_run.m pass 1 (band +-5 Hz, 400 m box, cap 1e3) -> cache/inv300.mat.
%  pass 2, per frequency: read Cs, form R0 and the methods of fe2_cfg (UQW1 form, gamma 1, s = gain_s(70))
%    on the stage-4 m grid, shots in chunks of 31 as fld_run.m -> cache/images_fe2.mat.
% mode 'smoke': only 15 Hz, tag _smoke; pass 2 uses the stage-4 1/snr for both noise slots (timing only,
% the smoke image is not evaluated). Memory gate before every LU (fe2_cfg mem_need_gb).
if nargin<1,mode='full';end
smoke=strcmp(mode,'smoke');[C,P]=fe2_cfg();maxNumCompThreads(C.threads);
ns=P.nshot;N=P.nz*P.nx;nf=numel(P.f);[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([P.nz P.nx],Zi(:),Xi(:));nin=numel(lin);
irec=sub2ind([P.nz P.nx],P.rec_iz,P.rec_ix);isrc=sub2ind([P.nz P.nx],P.src_iz,P.src_ix);col=sub2ind([N ns],irec,P.shot_of);W=P.W;
[EPS,nper,nb]=flip_signs(C,P);
fprintf('block flip B = %d: traces per shot %d-%d, blocks per shot %d, signs +%d/-%d over %d shots\n',C.B,min(nper),max(nper),nb,nnz(EPS>0),nnz(EPS<0),ns);
if smoke,jl=find(abs(P.f-15)<1e-9);tag='_smoke';else,jl=1:nf;tag='';end
csd=fullfile(C.cache,['cs' tag]);if ~isfolder(csd),mkdir(csd);end
% ---------------- pass 1 ----------------
f1=fullfile(C.cache,['inv300' tag '.mat']);ck=fullfile(C.cache,['pass1' tag '.mat']);
if ~isfile(f1)
 I=zeros(nin,nf,'single');If=I;done=false(1,nf);nsv=0;
 if isfile(ck),S=load(ck);I=S.I;If=S.If;done=S.done;end
 S4=load(fullfile(C.s4cache,'state1.mat'),'I');I4=S4.I;clear S4                % stage-4 pass-1 images (check)
 for j=jl
  if done(j),continue;end
  memgate(C);[pb0,~]=memstat();t0=tic;f=P.f(j);[L,U,Pp,Qq,R]=lu(helm_op(P.vmig,f,P));
  b=zeros(N,ns);b(sub2ind([N ns],isrc,(1:ns)'))=-W(j)/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
  b=zeros(N,ns);b(col)=P.D(:,j);V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;Cs=Us.*single(V(lin,:));clear V
  I(:,j)=real(sum(Cs,2));wcs(csd,j,Cs);clear Cs
  b=zeros(N,ns);b(col)=P.D(:,j).*EPS;V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;If(:,j)=real(sum(Us.*single(V(lin,:)),2));
  clear L U Pp Qq R b V Us
  done(j)=true;nsv=nsv+1;if mod(nsv,5)==0||j==jl(end),save(ck,'I','If','done','-v7.3');end
  [pb1,pk]=memstat();rd=norm(I(:,j)-I4(:,j))/max(norm(I4(:,j)),realmin);
  fprintf('fe2 pass 1 f %5.2f Hz (%d/%d): %.0f s, private %.2f GB, peak private %.2f GB (before LU %.2f), rel diff to stage-4 image %.2g  %s\n',f,nnz(done(jl)),numel(jl),toc(t0),pb1,pk,pb0,rd,datestr(now,'HH:MM:SS'));
 end
 if smoke
  fprintf('smoke: f %.2f Hz, source weight W = %.3g; |I| max %.3g, stage-4 |I| max %.3g\n',P.f(jl),W(jl),max(abs(I(:,jl))),max(abs(I4(:,jl))));
 elseif all(done)
  bw=round(400/P.dx);iv=zeros(nin,nf,'single');
  for j=1:nf
   bnd=abs(P.f-P.f(j))<=5+1e-9;p=box(sum(I(:,bnd),2),P,bw);pf=box(sum(If(:,bnd),2),P,bw);
   iv(:,j)=min(pf./max(p-pf,realmin('single')),1e3);
  end
  R0=sum(I,2);S4=load(fullfile(C.s4cache,'inv.mat'),'R0');
  fprintf('R0 vs stage-4 R0: max|diff|/max|R0| = %.3g\n',max(abs(R0-S4.R0))/max(abs(S4.R0)));
  save(f1,'iv','R0','-v7.3');fprintf('inv300 saved %s\n',datestr(now));clear I If I4 iv
 else
  error('fe2_run: pass 1 incomplete');
 end
end
% ---------------- pass 2 ----------------
fo=fullfile(C.cache,['images_fe2' tag '.mat']);if isfile(fo),return;end
S1=load(fullfile(C.s4cache,'inv.mat'),'iv');IV.B1=S1.iv;clear S1
if smoke,IV.B300=IV.B1;else,S=load(f1,'iv');IV.B300=S.iv;clear S;end
TS=fe2_tstar(C,P,lin);ms=C.ms;nm=size(C.meth,1);
xmax=1.05*double(max(TS.abs(:)))/min(ms);xg=single(linspace(0,xmax,1201));dxg=double(xg(2));   % stage-4 grid
if max(TS.dlt(:))>max(TS.abs(:)),error('Delta t* above the table range');end
s70=gain_s(70);ck2=fullfile(C.cache,['pass2' tag '.mat']);Img=struct();done2=false(1,nf);
if isfile(ck2),S=load(ck2);Img=S.Img;done2=S.done2;end
fam={'u100','abs';'u25','abs';'va','dlt'};nsv=0;
for j=jl
 if done2(j),continue;end
 t0=tic;f=P.f(j);Cs=rcs(csd,j,nin,ns);
 T.u100=tb(uq_tab(f,xg,C.sig,P.fref));T.u25=tb(uq_tab(f,xg,C.sig,25));T.va=tb(va_tab(f,xg,C.sig,C.fv));
 acc=struct('R0',zeros(nin,1,'single'));for q=1:nm,acc.(C.meth{q,1})=zeros(nin,numel(ms),'single');end
 for c=1:ceil(ns/31)
  sh=(c-1)*31+1:min(c*31,ns);Cc=Cs(:,sh);acc.R0=acc.R0+real(sum(Cc,2));
  for k=1:numel(ms)
   key='';
   for a=1:size(fam,1)
    if ~strcmp(key,fam{a,2})                                           % interpolation weights once per t* key
     key=fam{a,2};u=TS.(key)(:,sh)/ms(k)/dxg;i0=min(floor(u),numel(xg)-2);ww=u-single(i0);i0=i0+1;
    end
    it=@(t) t(i0).*(1-ww)+t(i0+1).*ww;Tb=T.(fam{a,1});A1=it(Tb.A);C1=it(Tb.C);P1=it(Tb.P);
    for q=find(strcmp(C.meth(:,3),fam{a,1}))'
     iv=IV.(C.meth{q,4})(:,j);H=(A1+s70*P1.*C1+P1.*iv)./(1+s70*C1+iv);
     acc.(C.meth{q,1})(:,k)=acc.(C.meth{q,1})(:,k)+real(sum(H.*Cc,2));
    end
   end
  end
 end
 for q=fieldnames(acc)',if isfield(Img,q{1}),Img.(q{1})=Img.(q{1})+acc.(q{1});else,Img.(q{1})=acc.(q{1});end,end
 done2(j)=true;nsv=nsv+1;if mod(nsv,5)==0||j==jl(end),save(ck2,'Img','done2','-v7.3');end
 fprintf('fe2 pass 2 f %5.2f Hz (%d/%d): %.0f s  %s\n',f,nnz(done2(jl)),numel(jl),toc(t0),datestr(now,'HH:MM:SS'));
end
if ~all(done2(jl)),error('fe2_run: pass 2 incomplete');end
% image file for fld_eval.m: stage-4 fields (unchanged copies, so that fld_eval runs as in stage 4 and
% reproduces stage4_field.json) plus the new methods; R0 is this run's (checked against stage 4 above)
S4=load(fullfile(C.s4cache,'images.mat'));X=struct('R0',Img.R0);m4=S4.meta.methods;
for q=m4,X.(q{1})=S4.(q{1});end
for q=1:nm,X.(C.meth{q,1})=Img.(C.meth{q,1});end
X.meta=struct('methods',{[m4(:)' C.meth(:,1)']},'m',ms,'fcut',S4.meta.fcut,'new',{C.meth(:,1)'},'noise',{C.meth(:,4)'});
save(fo,'-struct','X','-v7.3');fprintf('images saved %s (%s)\n',fo,datestr(now));
end

function [EPS,nper,nb]=flip_signs(C,P)
% est_block.m rule on field geometry: within a shot, receivers in order of position, B consecutive share
% one sign; rng(flip_seed + B,'twister'), signs s(nb x nshot).
ns=P.nshot;nper=accumarray(P.shot_of(:),1,[ns 1]);nb=ceil(max(nper)/C.B);
rng(C.flip_seed+C.B,'twister');s=single(2*(rand(nb,ns)>.5)-1);blk=zeros(P.ntrace,1);
for k=1:ns
 t=find(P.shot_of==k);[~,o]=sort(P.rec_ix(t));r=zeros(numel(t),1);r(o)=1:numel(t);blk(t)=ceil(r/C.B);
end
EPS=s(sub2ind([nb ns],blk,P.shot_of(:)));EPS=EPS(:);
end

function T=tb(T),T.A=single(T.A);T.C=single(T.C);T.P=single(T.P);end

function wcs(d,j,Cs)
fid=fopen(fullfile(d,sprintf('f%03d.bin',j)),'w');fwrite(fid,[real(Cs(:)) imag(Cs(:))].','single');fclose(fid);
end

function Cs=rcs(d,j,nin,ns)
fid=fopen(fullfile(d,sprintf('f%03d.bin',j)));a=fread(fid,[2 nin*ns],'*single');fclose(fid);
Cs=reshape(complex(a(1,:),a(2,:)),nin,ns);
end

function y=box(v,P,bw)
y=reshape(v,P.nzi,P.nxi).^2;y=movmean(movmean(y,bw,1),bw,2);y=y(:);
end

function [cur,pk]=memstat()
% private bytes of this process now and at peak (GB), via .NET
cur=NaN;pk=NaN;
try,p=System.Diagnostics.Process.GetCurrentProcess();p.Refresh();cur=double(p.PrivateMemorySize64)/2^30;pk=double(p.PeakPagedMemorySize64)/2^30;catch,end
end

function memgate(C)
% wait until the system has mem_need_gb of commit headroom before an LU (other users' processes)
if isnan(C.mem_need_gb),return;end
t0=tic;last=-Inf;
while true
 [~,sv]=memory;av=sv.SystemMemory.Available/2^30;
 if av>=C.mem_need_gb,break;end
 if toc(t0)-last>=600,fprintf('memgate: %.2f GB available < %.2f GB needed, waiting (%s)\n',av,C.mem_need_gb,datestr(now,'HH:MM:SS'));last=toc(t0);end
 pause(60);
end
end
