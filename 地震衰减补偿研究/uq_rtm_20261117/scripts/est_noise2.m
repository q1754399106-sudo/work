function est_noise2(C)
% Pass 1 for stage 2 (as est_noise.m, round 3): band-summed random-polarity 1/snr(x,f) for every
% dataset of C.sets, migration velocity M.vmig*C.vscale. Output cache/inv_<name><tag>.mat.
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model(C.name,P);vmig=M.vmig*C.vscale;
fo=fullfile(P.cache,['inv_' C.name C.tag '.mat']);if isfile(fo),return;end
D=load(fullfile(P.cache,['data_' C.name '.mat']));ns=numel(P.ixs);nz=P.nz;nx=P.nx;N=nz*nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([nz nx],Zi(:),Xi(:));nin=numel(lin);
irec=sub2ind([nz nx],P.izs*ones(size(P.ixr)),P.ixr);isrc=sub2ind([nz nx],P.izs*ones(size(P.ixs)),P.ixs);
rng(P.flip_seed,'twister');EPS=single(2*(rand(numel(P.ixr),ns)>.5)-1);
t=(1:P.fd.nt)'*P.fd.dt;a=(pi*P.fd.fpeak*(t-P.fd.t0)).^2;q=(1-2*a).*exp(-a);
sets=C.sets(:,1);out=fullfile(P.cache,['est_' C.name C.tag]);if ~isfolder(out),mkdir(out);end
for j=1:numel(P.f)
 fn=fullfile(out,sprintf('f%03d.mat',j));if isfile(fn),continue;end
 t0=tic;f=P.f(j);W=(1i*2*pi*f)*(q.'*exp(-1i*2*pi*t*f))*P.fd.dt;
 A=helm_op(vmig,f,P);[L,U,Pp,Qq,R]=lu(A);clear A
 b=zeros(N,ns);b(sub2ind([N ns],isrc,1:ns))=-W/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
 S=struct();
 for d=1:numel(sets)
  Dd=squeeze(D.(sets{d})(:,j,:));
  b=zeros(N,ns);b(irec,:)=Dd;V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;S.(sets{d})=real(sum(Us.*single(V(lin,:)),2));
  b(irec,:)=Dd.*EPS;V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;S.([sets{d} '_f'])=real(sum(Us.*single(V(lin,:)),2));
 end
 save(fn,'-struct','S');fprintf('est %s%s f %4.1f Hz (%d/%d): %.0f s\n',C.name,C.tag,f,j,numel(P.f),toc(t0));
end
bw=round(P.snr_box/P.dx);INV=struct();
for d=1:numel(sets)
 I=zeros(nin,numel(P.f),'single');If=I;
 for j=1:numel(P.f)
  S=load(fullfile(out,sprintf('f%03d.mat',j)),sets{d},[sets{d} '_f']);I(:,j)=S.(sets{d});If(:,j)=S.([sets{d} '_f']);
 end
 iv=zeros(nin,numel(P.f),'single');
 for j=1:numel(P.f)
  bnd=abs(P.f-P.f(j))<=P.snr_band+1e-9;
  p=box(sum(I(:,bnd),2),P,bw);pf=box(sum(If(:,bnd),2),P,bw);
  iv(:,j)=min(pf./max(p-pf,realmin('single')),1/P.snr_min);
 end
 INV.(sets{d})=iv;fprintf('inv %s%s %s: median 1/snr at 10/25/40 Hz: %.3g %.3g %.3g\n',C.name,C.tag,sets{d}, ...
  median(iv(:,P.f==10)),median(iv(:,P.f==25)),median(iv(:,P.f==40)));
end
save(fo,'-struct','INV','-v7.3');
end

function y=box(v,P,bw)
y=reshape(v,P.nzi,P.nxi).^2;y=movmean(movmean(y,bw,1),bw,2);y=y(:);
end
