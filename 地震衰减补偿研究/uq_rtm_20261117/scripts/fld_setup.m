function P=fld_setup()
% Stage-4 constants and data (PLAN_STAGE4_ZH.md), read from cache/fld (fld_prep.py).
P.root=fileparts(fileparts(mfilename('fullpath')));P.cache=fullfile(P.root,'cache','fld');
M=jsondecode(fileread(fullfile(P.cache,'meta.json')));
P.dx=M.dx;P.npml=M.npml;P.nzi=M.nzi;P.nxi=M.nxi;P.nz=M.nz;P.nx=M.nx;
P.x=M.xmin+((1:P.nx)-P.npml-1)*P.dx;P.z=((1:P.nz)-P.npml-1)*P.dx;P.in={P.npml+(1:P.nzi),P.npml+(1:P.nxi)};
P.f=M.f(:)';P.nshot=M.nshot;P.ntrace=M.ntrace;P.threads=12;P.fref=100;P.sig=.447;
P.shot_of=double(rdc(P,'shot_of.bin','int32'));P.rec_ix=double(rdc(P,'rec_ix.bin','int32'));P.rec_iz=double(rdc(P,'rec_iz.bin','int32'));
P.src_ix=double(rdc(P,'src_ix.bin','int32'));P.src_iz=double(rdc(P,'src_iz.bin','int32'));
d=rdc(P,'D.bin','single');d=reshape(d,2,numel(P.f),P.ntrace);P.D=single(squeeze(complex(d(1,:,:),d(2,:,:))).');   % ntrace x nf
v=rdc(P,'vp.bin','single');P.vp=double(reshape(v,P.nx,P.nz).');
P.vmig=gsm(P.vp,50/P.dx);P.vmig(P.vp<1500)=P.vp(P.vp<1500);                % keep water velocity exact
P.Qm=120*ones(P.nz,P.nx);P.Qm(P.vp<1500)=Inf;                               % Q model: water lossless
f=P.f;w=ones(size(f));a=f<6;w(a)=.5*(1-cos(pi*(f(a)-4)/2));a=f>23;w(a)=.5*(1+cos(pi*(f(a)-23)/2));P.W=w;
P.src_x=P.x(P.src_ix);
end

function v=rdc(P,n,t)
fid=fopen(fullfile(P.cache,n));v=fread(fid,Inf,['*' t]);fclose(fid);
end

function y=gsm(x,s)
h=exp(-.5*((-ceil(3*s):ceil(3*s))/s).^2);h=h/sum(h);m=(numel(h)-1)/2;
xp=[repmat(x(1,:),m,1);x;repmat(x(end,:),m,1)];xp=[repmat(xp(:,1),1,m) xp repmat(xp(:,end),1,m)];
y=conv2(h,h,xp,'valid');
end
