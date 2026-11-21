function va_mig(name,jlist)
% Step 1 of the VA-UQ document (PLAN_ZH.md): migrate the m30 data of a development model (D2, D4)
% with the in-band phase velocity v_band (25 Hz), and form R0 and the methods of va_cfg on the E1 grid.
% jlist (optional): frequency indices only (smoke test, tag _smoke). Resumable: state every 5 f.
% Output cache/images_va_<name><tag>.mat.
[C,P]=va_cfg(name);if nargin<2||isempty(jlist),jlist=1:numel(P.f);else,C.tag='_smoke';end
maxNumCompThreads(P.threads);
M=make_model(name,P);[vb,vchk]=va_vband(M,P,C);
fprintf('v_band %s: smoothing check max|gsmooth(v)-vmig| = %.3g m/s; vb/vmig in [%.4f, %.4f]\n',name,vchk,min(vb(:)./M.vmig(:)),max(vb(:)./M.vmig(:)));
D=load(fullfile(C.src_cache,['data_' name '.mat']),C.dataset);D=D.(C.dataset);
INV=load(fullfile(C.src_cache,['inv_' name '.mat']),C.dataset);INV=INV.(C.dataset);   % est_noise2 B = 1, v100 (PLAN 4)
ns=numel(P.ixs);nz=P.nz;nx=P.nx;N=nz*nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([nz nx],Zi(:),Xi(:));
irec=sub2ind([nz nx],P.izs*ones(size(P.ixr)),P.ixr);isrc=sub2ind([nz nx],P.izs*ones(size(P.ixs)),P.ixs);
TS=va_tstar(M,P,C,vb,lin,Xi(:));
xmax=1.05*double(max(TS.abs(:)))/min(C.m);xg=single(linspace(0,xmax,1201));dxg=double(xg(2));
t=(1:P.fd.nt)'*P.fd.dt;a=(pi*P.fd.fpeak*(t-P.fd.t0)).^2;q=(1-2*a).*exp(-a);
out=fullfile(P.cache,['vamig_' name C.tag]);if ~isfolder(out),mkdir(out);end
ck=fullfile(out,'state.mat');I=struct();done=[];
if isfile(ck),S=load(ck);I=S.I;done=S.done;end
s70=gain_s(70);nsv=0;
for j=jlist
 if ismember(j,done),continue;end
 t0=tic;f=P.f(j);w=2*pi*f;W=(1i*w)*(q.'*exp(-1i*w*t))*P.fd.dt;
 A=helm_op(vb,f,P);[L,U,Pp,Qq,R]=lu(A);clear A
 b=zeros(N,ns);b(sub2ind([N ns],isrc,1:ns))=-W/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*w^2;
 Dd=squeeze(D(:,j,:));b=zeros(N,ns);b(irec,:)=Dd;V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;
 Cs=Us.*single(V(lin,:));clear V b Us L U Pp Qq R
 img=struct('R0',real(sum(Cs,2)));iv=INV(:,j);
 tab=va_tables(f,xg,C,P);
 for r=1:numel(C.m)
  for mm=C.methods
   T=tab.(mm{1});u=TS.(T.key)/C.m(r)/dxg;i0=min(floor(u),numel(xg)-2);ww=u-single(i0);i0=i0+1;
   it=@(tb) tb(i0).*(1-ww)+tb(i0+1).*ww;
   if T.phase_only,H=it(T.P);
   else,A1=it(T.A);C1=it(T.C);P1=it(T.P);H=(A1+s70*P1.*C1+P1.*iv)./(1+s70*C1+iv);end    % UQW1 form, gamma 1
   img.(mm{1})(:,r)=real(sum(H.*Cs,2));
  end
 end
 clear Cs H A1 C1 P1 u i0 ww
 I=accum(I,img);done(end+1)=j;nsv=nsv+1; %#ok<AGROW>
 if mod(nsv,5)==0||j==jlist(end),save(ck,'I','done','-v7.3');end
 mu=NaN;try,[uu,~]=memory;mu=uu.MemUsedMATLAB/2^30;catch,end
 fprintf('vamig %s%s f %4.1f Hz (%d/%d): %.0f s, MATLAB mem %.2f GB  %s\n',name,C.tag,f,numel(done),numel(jlist),toc(t0),mu,datestr(now,'HH:MM:SS'));
end
save(ck,'I','done','-v7.3');
if all(ismember(jlist,done))
 I.meta=struct('m',C.m,'methods',{C.methods},'freq',P.f(jlist),'cfg',C,'xmax',xmax);
 save(fullfile(P.cache,['images_va_' name C.tag '.mat']),'-struct','I','-v7.3');
 fprintf('vamig %s%s done %s\n',name,C.tag,datestr(now));
end
end

function tab=va_tables(f,xg,C,P)
% Expectation tables per method: key = t* field (abs: absolute two-way t*, dlt: anchored Delta t*)
fr=P.fd.fref;                                                          % 100 Hz (GSLS reference)
tab.PH100=tt(uq_tab(f,xg,0,fr),'abs',true);
tab.UQW1_100=tt(uq_tab(f,xg,C.sig,fr),'abs',false);
tab.UQW1_25=tt(uq_tab(f,xg,C.sig,C.f_band),'abs',false);
tab.VA=tt(va_tab(f,xg,C.sig,C.fv),'dlt',false);
tab.PH25=tt(uq_tab(f,xg,0,C.f_band),'abs',true);
tab.VAfix=tt(va_tab(f,xg,C.sig,C.f_band),'dlt',false);
end

function T=tt(T,key,po)
T.A=single(T.A);T.C=single(T.C);T.P=single(T.P);T.key=key;T.phase_only=po;
end

function I=accum(I,S)
for q=fieldnames(S)'
 if ~isfield(I,q{1}),I.(q{1})=S.(q{1});else,I.(q{1})=I.(q{1})+S.(q{1});end
end
end
