function est_block(name,jlist,otag)
% Block-randomised polarity noise estimate (PLAN_BLOCKFLIP_ZH.md, development on D1-D4).
% As est_noise2.m, but the random sign is shared by B consecutive receivers of a shot:
%   B = 1   identical signs to est_noise2 (rng(P.flip_seed)), reproduction check
%   B = 8, 32 (64 m, 256 m)      B = 300 one sign per shot
% Datasets m30, m40, g30 (data_<name>.mat) and l30 (laterally correlated noise generated here for the D
% models with dev seeds, data_<name>_l30.mat). Truth: image of the noise alone = I(X) - I(att).
% Per frequency stores unflipped I(X), I(att) and flipped images; then band sums (+-P.snr_band) and the
% 400 m box exactly as est_noise2. Output cache/invblk_<name><otag>.mat with fields <set>_B<B> and
% <set>_true (1/snr, nin x nf). jlist (optional) restricts frequencies (smoke test), otag output tag.
P=rtm_setup();maxNumCompThreads(P.threads);M=make_model(name,P);vmig=M.vmig;
if nargin<2||isempty(jlist),jlist=1:numel(P.f);end
if nargin<3,otag='';end
BS=[1 8 32 300];sets={'m30','m40','g30','l30'};
fo=fullfile(P.cache,['invblk_' name otag '.mat']);if isfile(fo),return;end
D=load(fullfile(P.cache,['data_' name '.mat']),'att','m30','m40','g30');
fl=fullfile(P.cache,['data_' name '_l30.mat']);if ~isfile(fl),make_l30(name,P,fl);end
L30=load(fl);D.l30=L30.l30;clear L30
ns=numel(P.ixs);nr=numel(P.ixr);nz=P.nz;nx=P.nx;N=nz*nx;
[Zi,Xi]=ndgrid(P.in{1},P.in{2});lin=sub2ind([nz nx],Zi(:),Xi(:));nin=numel(lin);
irec=sub2ind([nz nx],P.izs*ones(size(P.ixr)),P.ixr);isrc=sub2ind([nz nx],P.izs*ones(size(P.ixs)),P.ixs);
EPS=cell(1,numel(BS));
for b=1:numel(BS)
 B=BS(b);
 if B==1,rng(P.flip_seed,'twister');EPS{b}=single(2*(rand(nr,ns)>.5)-1);
 else,rng(P.flip_seed+B,'twister');nb=ceil(nr/B);s=single(2*(rand(nb,ns)>.5)-1);EPS{b}=s(ceil((1:nr)/B),:);end
end
t=(1:P.fd.nt)'*P.fd.dt;a=(pi*P.fd.fpeak*(t-P.fd.t0)).^2;q=(1-2*a).*exp(-a);
out=fullfile(P.cache,['estblk_' name otag]);if ~isfolder(out),mkdir(out);end
for j=jlist
 fn=fullfile(out,sprintf('f%03d.mat',j));if isfile(fn),continue;end
 t0=tic;f=P.f(j);W=(1i*2*pi*f)*(q.'*exp(-1i*2*pi*t*f))*P.fd.dt;
 A=helm_op(vmig,f,P);[L,U,Pp,Qq,R]=lu(A);clear A
 b=zeros(N,ns);b(sub2ind([N ns],isrc,1:ns))=-W/P.dx^2;Us=Qq*(U\(L\(Pp*(R\b))));Us=conj(single(Us(lin,:)))*(2*pi*f)^2;
 img=@(Dd) mig1(Dd,Us,L,U,Pp,Qq,R,N,ns,irec,lin,P);
 S=struct();S.att=img(squeeze(D.att(:,j,:)));
 for d=1:numel(sets)
  Dd=squeeze(D.(sets{d})(:,j,:));S.(sets{d})=img(Dd);
  for k=1:numel(BS),S.(sprintf('%s_B%d',sets{d},BS(k)))=img(Dd.*EPS{k});end
 end
 save(fn,'-struct','S');fprintf('estblk %s%s f %4.1f Hz (%d/%d): %.0f s\n',name,otag,f,j,numel(P.f),toc(t0));
end
if ~isequal(jlist,1:numel(P.f)),return;end                 % smoke: per-frequency files only
bw=round(P.snr_box/P.dx);nf=numel(P.f);INV=struct();
Ia=zeros(nin,nf,'single');for j=1:nf,S=load(fullfile(out,sprintf('f%03d.mat',j)),'att');Ia(:,j)=S.att;end
for d=1:numel(sets)
 sd=sets{d};I=zeros(nin,nf,'single');IB=zeros(nin,nf,numel(BS),'single');
 for j=1:nf
  S=load(fullfile(out,sprintf('f%03d.mat',j)));I(:,j)=S.(sd);
  for k=1:numel(BS),IB(:,j,k)=S.(sprintf('%s_B%d',sd,BS(k)));end
 end
 In=I-Ia;tv=zeros(nin,nf,'single');iv=zeros(nin,nf,numel(BS),'single');
 for j=1:nf
  bnd=abs(P.f-P.f(j))<=P.snr_band+1e-9;p=box(sum(I(:,bnd),2),P,bw);
  tv(:,j)=min(box(sum(In(:,bnd),2),P,bw)./max(box(sum(Ia(:,bnd),2),P,bw),realmin('single')),1/P.snr_min);
  for k=1:numel(BS)
   pf=box(sum(IB(:,bnd,k),2),P,bw);iv(:,j,k)=min(pf./max(p-pf,realmin('single')),1/P.snr_min);
  end
 end
 INV.([sd '_true'])=tv;for k=1:numel(BS),INV.(sprintf('%s_B%d',sd,BS(k)))=iv(:,:,k);end
end
INV.BS=BS;save(fo,'-struct','INV','-v7.3');fprintf('invblk %s%s saved\n',name,otag);
end

function y=mig1(Dd,Us,L,U,Pp,Qq,R,N,ns,irec,lin,P)
b=zeros(N,ns);b(irec,:)=Dd;V=-(R\(Pp'*(L'\(U'\(Qq'*b)))))/P.dx^2;y=real(sum(Us.*single(V(lin,:)),2));
end

function y=box(v,P,bw)
y=reshape(v,P.nzi,P.nxi).^2;y=movmean(movmean(y,bw,1),bw,2);y=y(:);
end

function make_l30(name,P,fl)
% Laterally correlated noise for a D model, recipe of prep_data.m (stage 3), seed P.lc.seed + 1000 + model no.
TD=load(fullfile(P.cache,['tdata_' name '.mat']),'att');ATT=TD.att;[nt,nr,ns]=size(ATT);
off=(P.ixr(:)-P.ixs(:)')*P.dx;tw=P.t>=1.2&P.t<=1.8;near=abs(off)<=600;
a=reshape(ATT(tw,:,:),nnz(tw),[]);rms_s=sqrt(mean(a(:,near(:)).^2,'all'));
fq=(0:nt-1)'/(nt*P.dt);fq=min(fq,1/P.dt-fq);band=ones(nt,1);
band(fq<3|fq>55)=0;b=fq>=3&fq<6;band(b)=.5*(1-cos(pi*(fq(b)-3)/3));b=fq>50&fq<=55;band(b)=.5*(1+cos(pi*(fq(b)-50)/5));
rng(P.lc.seed+1000+str2double(name(2:end)),'twister');n=randn(nt,nr*ns);n=real(ifft(fft(n).*band));n=reshape(n,nt,nr,ns);
h=exp(-.5*((-3*P.lc.sigma_tr:3*P.lc.sigma_tr)/P.lc.sigma_tr).^2);h=h/sum(h);
for k=1:ns,n(:,:,k)=conv2(1,h,n(:,:,k),'same');end
n=single(n/sqrt(mean(n(:).^2))*rms_s/10^(P.lc.db/20));X=ATT+n;
f=3:.5:55;E=single(exp(-1i*2*pi*P.t*f));Df=reshape(E.'*reshape(X,nt,nr*ns),[],nr,ns);l30=permute(Df,[2 1 3]); %#ok<NASGU>
save(fl,'l30','-v7.3');fprintf('make_l30 %s: signal rms %.3g, seed %d\n',name,rms_s,P.lc.seed+1000+str2double(name(2:end)));
end
