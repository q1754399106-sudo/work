function migrate(name)
% Frequency-domain two-way (reverse-time) migration with per-point compensation operators in the
% imaging condition (route A), PLAN_STAGE1_ZH.md. For each frequency: one LU of the acoustic
% operator in the migration velocity; source fields U_s; adjoint (back-propagated) receiver fields
% V_s of every dataset; per-shot image contribution C_s = w^2 conj(U_s) V_s on the interior. Images
%   I = sum_f sum_s Re[ H(x_hat(s,x), 1/snr(x,f), f) C_s(x,f) ],
% x_hat = two-way t* of the compensation Q model (two_way.m) divided by the Q multiplier m.
% 1/snr(x,f): image-domain, from the uncompensated even/odd-shot partial images of the same
% dataset at the same frequency, snr = 4 <Ie Io> / <(Ie-Io)^2> in a P.snr_box square window.
% Output cache/images_<name>.mat; resumable (running sums saved every 10 frequencies, cache/mig_<name>/state.mat).
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model(name,P);F=tstar_fields(M,P);
D=load(fullfile(P.cache,['data_' name '.mat']));ns=numel(P.ixs);nz=P.nz;nx=P.nx;N=nz*nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([nz nx],Zi(:),Xi(:));nin=numel(lin);
TS=zeros(nin,ns,'single');
for k=1:ns,[~,tt]=two_way(F,k,P);TS(:,k)=tt(lin);end
clear F
irec=sub2ind([nz nx],P.izs*ones(size(P.ixr)),P.ixr);isrc=sub2ind([nz nx],P.izs*ones(size(P.ixs)),P.ixs);
W=wsrc(P);out=fullfile(P.cache,['mig_' name]);if ~isfolder(out),mkdir(out);end
% methods: name, sig (0 or P.sig), gamma, uses 1/snr, gain (dB) of the floor s
meth={'PH',0,NaN,0,NaN;'IQ20',0,1,0,20;'IQ30',0,1,0,30;'IQ40',0,1,0,40;'IQ50',0,1,0,50;'IQ70',0,1,0,70; ...
 'W0',0,0,1,70;'UQW0',P.sig,0,1,70;'UQW1',P.sig,1,1,70;'UQIQ',P.sig,1,0,70};
full={'att',[.7 1 1.4];'n00',P.m;'m12',P.m};                % datasets with operators, and their m
xmax=1.05*double(max(TS(:)))/min(P.m);xg=single(linspace(0,xmax,1201));dxg=double(xg(2));
bw=round(P.snr_box/P.dx);
cdata=fieldnames(D);cdata=cdata(startsWith(cdata,'C_'));
ck=fullfile(out,'state.mat');I=struct();jdone=0;
if isfile(ck),S=load(ck);I=S.I;jdone=S.jdone;end
for j=jdone+1:numel(P.f)
 t0=tic;f=P.f(j);w=2*pi*f;
 A=helm_op(M.vmig,f,P);[L,U,Pp,Qq,R]=lu(A);clear A
 solve=@(b) Qq*(U\(L\(Pp*(R\b))));solveH=@(b) R\(Pp'*(L'\(U'\(Qq'*b))));
 b=zeros(N,ns);b(sub2ind([N ns],isrc,1:ns))=-W(j)/P.dx^2;Us=solve(b);Us=conj(single(Us(lin,:)))*w^2;
 img=struct();
 sets=[{'ref'};full(:,1);cdata];
 for d=1:numel(sets)
  Dd=squeeze(D.(sets{d})(:,j,:));                           % nrec x ns
  b=zeros(N,ns);b(irec,:)=Dd;V=-solveH(b)/P.dx^2;C=Us.*single(V(lin,:));clear V b
  img.(sets{d}).R0=real(sum(C,2));
  k=find(strcmp(full(:,1),sets{d}));if isempty(k),continue;end
  Ie=reshape(real(sum(C(:,1:2:end),2)),P.nzi,P.nxi);Io=reshape(real(sum(C(:,2:2:end),2)),P.nzi,P.nxi);
  sn=4*max(box(Ie.*Io,bw),0)./max(box((Ie-Io).^2,bw),realmin('single'));inv=single(1./max(sn(:),P.snr_min));
  img.(sets{d}).inv_med=median(inv);
  T0=uq_tab(f,xg,0,P.fd.fref);T1=uq_tab(f,xg,P.sig,P.fd.fref);
  for mi=1:numel(full{k,2})
   m=full{k,2}(mi);u=TS/m/dxg;i0=min(floor(u),numel(xg)-2);ww=u-single(i0);i0=i0+1;
   it=@(t) t(i0).*(1-ww)+t(i0+1).*ww;
   A0=it(single(T0.A));C0=it(single(T0.C));P0=it(single(T0.P));
   A1=it(single(T1.A));C1=it(single(T1.C));P1=it(single(T1.P));
   for q=1:size(meth,1)
    [nm,sg,ga,us,g]=meth{q,:};
    if sg==0,Aa=A0;Cc=C0;Pq=P0;else,Aa=A1;Cc=C1;Pq=P1;end
    if strcmp(nm,'PH'),H=Pq;
    else
     s=gain_s(g);iv=0;if us,iv=inv;end
     H=(Aa+ga*(s*Pq.*Cc+Pq.*iv))./(1+s*Cc+iv);
    end
    img.(sets{d}).(nm)(:,mi)=real(sum(H.*C,2));
   end
  end
 end
 I=accum(I,img,j);
 if mod(j,10)==0||j==numel(P.f),jdone=j;save(ck,'I','jdone','-v7.3');end
 fprintf('%s f %4.1f Hz (%d/%d): %.0f s\n',name,f,j,numel(P.f),toc(t0));
end
I.meta=struct('methods',{meth(:,1)'},'m_full',{full},'nzi',P.nzi,'nxi',P.nxi,'f',P.f);
save(fullfile(P.cache,['images_' name '.mat']),'-struct','I','-v7.3');
end

function I=accum(I,S,j)
for d=fieldnames(S)'
 for q=fieldnames(S.(d{1}))'
  if strcmp(q{1},'inv_med'),I.(d{1}).inv_med(j)=S.(d{1}).inv_med;continue;end
  if ~isfield(I,d{1})||~isfield(I.(d{1}),q{1}),I.(d{1}).(q{1})=S.(d{1}).(q{1});else,I.(d{1}).(q{1})=I.(d{1}).(q{1})+S.(d{1}).(q{1});end
 end
end
end

function y=box(x,b)
y=movmean(movmean(x,b,1),b,2);
end

function W=wsrc(P)
t=(1:P.fd.nt)'*P.fd.dt;a=(pi*P.fd.fpeak*(t-P.fd.t0)).^2;q=(1-2*a).*exp(-a);
W=(1i*2*pi*P.f).*(q.'*exp(-1i*2*pi*t*P.f))*P.fd.dt;
end
