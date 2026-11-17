function migrate2(C)
% Stage-2 migration (PLAN_STAGE2_ZH.md): route A with the band-summed noise estimate of est_noise2.
% For each dataset of C.sets and each (Q-error form, m) combination, images of C.methods; plus R0
% of every dataset and of ref. Output cache/images2_<name><tag>.mat; resumable (state every 10 f).
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model(C.name,P);vmig=M.vmig*C.vscale;
if ~isfield(C,'itag'),C.itag=C.tag;end
if ~isfield(C,'fixed_methods'),C2=s2_cfg(C.name,C.vscale);C.methods=C2.methods;end   % addendum: methods read at start
D=load(fullfile(P.cache,['data_' C.name '.mat']));INV=load(fullfile(P.cache,['inv_' C.name C.itag '.mat']));
ns=numel(P.ixs);nz=P.nz;nx=P.nx;N=nz*nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([nz nx],Zi(:),Xi(:));
irec=sub2ind([nz nx],P.izs*ones(size(P.ixr)),P.ixr);isrc=sub2ind([nz nx],P.izs*ones(size(P.ixs)),P.ixs);
% combinations: per dataset, rows {form label, m, t* key, scale}
keys={};TS=struct();combo=struct();
for d=1:size(C.sets,1)
 rows={};F=C.sets{d,2};
 for r=1:size(F,1)
  for m=F{r,2}
   if strcmp(F{r,1},'E4'),key=sprintf('E4_%03d',round(100*m));sc=1;else,key=F{r,1};sc=1/m;end
   rows(end+1,:)={F{r,1},m,key,sc}; %#ok<AGROW>
   if ~ismember(key,keys),keys{end+1}=key;end %#ok<AGROW>
  end
 end
 combo.(C.sets{d,1})=rows;
end
for k=1:numel(keys)
 key=keys{k};form=key(1:2);mq=1;if strcmp(form,'E4'),mq=str2double(key(4:end))/100;end
 if strcmp(key,'E1'),tg=C.itag;else,tg=[C.itag '_' key];end
 F=tstar_fields(M,P,make_qc(M,P,form,mq),vmig,tg);T=zeros(numel(lin),ns,'single');
 for s=1:ns,[~,tt]=two_way(F,s,P);T(:,s)=tt(lin);end
 TS.(key)=T;clear F
end
xmax=0;for k=1:numel(keys),xmax=max(xmax,double(max(TS.(keys{k})(:))));end
xmax=1.05*xmax/min(P.m);xg=single(linspace(0,xmax,1201));dxg=double(xg(2));
% LP cut-off per dataset: highest f of the run from 10 Hz where the deep-window median snr >= 1
z=P.z(P.in{1});x=P.x(P.in{2});wd=(z>=820&z<=1550)'&(x>=300&x<=2100);fc=struct();
for d=1:size(C.sets,1)
 iv=INV.(C.sets{d,1});sm=1./median(iv(wd(:),:),1);j0=find(P.f>=10,1);j=j0;
 while j<numel(P.f)&&sm(j+1)>=1,j=j+1;end
 if sm(j0)<1,j=j0;end
 fc.(C.sets{d,1})=P.f(j);fprintf('LP cut-off %s%s %s: %.1f Hz\n',C.name,C.tag,C.sets{d,1},P.f(j));
end
t=(1:P.fd.nt)'*P.fd.dt;a=(pi*P.fd.fpeak*(t-P.fd.t0)).^2;q=(1-2*a).*exp(-a);
out=fullfile(P.cache,['mig2_' C.name C.tag]);if ~isfolder(out),mkdir(out);end
ck=fullfile(out,'state.mat');I=struct();jdone=0;
if isfile(ck),S=load(ck);I=S.I;jdone=S.jdone;end
s30=gain_s(30);s70=gain_s(70);
for j=jdone+1:numel(P.f)
 t0=tic;f=P.f(j);w=2*pi*f;W=(1i*w)*(q.'*exp(-1i*w*t))*P.fd.dt;
 A=helm_op(vmig,f,P);[L,U,Pp,Qq,R]=lu(A);clear A
 b=zeros(N,ns);b(sub2ind([N ns],isrc,1:ns))=-W/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*w^2;
 T0=uq_tab(f,xg,0,P.fd.fref);T1=uq_tab(f,xg,P.sig,P.fd.fref);
 img=struct();
 for d=[{'ref'};C.sets(:,1)]'
  Dd=squeeze(D.(d{1})(:,j,:));b=zeros(N,ns);b(irec,:)=Dd;V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;
  Cs=Us.*single(V(lin,:));clear V b
  img.(d{1}).R0=real(sum(Cs,2));if strcmp(d{1},'ref'),continue;end
  iv=INV.(d{1})(:,j);fcut=fc.(d{1});Lf=single(f<=fcut-5)+single(f>fcut-5&f<fcut)*.5*(1+cos(pi*(f-fcut+5)/5));
  rows=combo.(d{1});
  for r=1:size(rows,1)
   u=TS.(rows{r,3})*rows{r,4}/dxg;i0=min(floor(u),numel(xg)-2);ww=u-single(i0);i0=i0+1;
   it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
   A0=it(single(T0.A));C0=it(single(T0.C));P0=it(single(T0.P));A1=it(single(T1.A));C1=it(single(T1.C));P1=it(single(T1.P));
   for mm=C.methods
    switch mm{1}
     case 'PH',H=P0;
     case 'IQ30',H=(A0+s30*P0.*C0)./(1+s30*C0);
     case 'IQ70',H=(A0+s70*P0.*C0)./(1+s70*C0);
     case 'LP',H=A0*Lf;
     case 'W0',H=A0./(1+s70*C0+iv);
     case 'UQW0',H=A1./(1+s70*C1+iv);
     case 'W1',H=(A0+s70*P0.*C0+P0.*iv)./(1+s70*C0+iv);
     case 'UQW1',H=(A1+s70*P1.*C1+P1.*iv)./(1+s70*C1+iv);
    end
    img.(d{1}).(mm{1})(:,r)=real(sum(H.*Cs,2));
   end
  end
 end
 I=accum(I,img);
 if mod(j,10)==0||j==numel(P.f),jdone=j;save(ck,'I','jdone','-v7.3');end
 fprintf('mig2 %s%s f %4.1f Hz (%d/%d): %.0f s\n',C.name,C.tag,f,j,numel(P.f),toc(t0));
end
I.meta=struct('combo',combo,'methods',{C.methods},'fc',fc,'tag',C.tag,'vscale',C.vscale);
save(fullfile(P.cache,['images2_' C.name C.tag '.mat']),'-struct','I','-v7.3');
end

function I=accum(I,S)
for d=fieldnames(S)'
 for q=fieldnames(S.(d{1}))'
  if ~isfield(I,d{1})||~isfield(I.(d{1}),q{1}),I.(d{1}).(q{1})=S.(d{1}).(q{1});else,I.(d{1}).(q{1})=I.(d{1}).(q{1})+S.(d{1}).(q{1});end
 end
end
end
