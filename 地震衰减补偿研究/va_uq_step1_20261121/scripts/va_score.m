function R=va_score(name,tag)
% Scores of PLAN_ZH.md section 5 for one model: S (E1 metric of score_stage2.m, identical code) of R0
% and every method on the m grid, against the reference = elastic data migrated with v100
% (uq_rtm_20261117 cache/images2_<name>.mat, field ref.R0). R0_old = m30 migrated with v100 (same file).
% Writes analysis/va_step1_<name><tag>.json.
if nargin<2,tag='';end
[C,P]=va_cfg(name);I=load(fullfile(P.cache,['images_va_' name tag '.mat']));
O=load(fullfile(C.src_cache,['images2_' name '.mat']),'ref','m30');
z=P.z(P.in{1});x=P.x(P.in{2});ref=prep(O.ref.R0,P);
R=struct('model',name,'tag',tag,'m',I.meta.m,'anchor_iface',C.anchor_iface, ...
 'anchor_depth_range_m',[min(anchor(C,name,P)) max(anchor(C,name,P))],'fv',C.fv,'n_freq',numel(I.meta.freq));
for wn=fieldnames(C.win)'
 W=C.win.(wn{1});T=struct();
 T.R0=score(prep(I.R0,P),ref,z,x,W);T.R0_old=score(prep(O.m30.R0,P),ref,z,x,W);
 for mm=I.meta.methods
  s=zeros(1,numel(I.meta.m));for r=1:numel(s),s(r)=score(prep(I.(mm{1})(:,r),P),ref,z,x,W);end
  T.(mm{1})=s;T.([mm{1} '_mean'])=mean(s);
 end
 R.(wn{1})=T;
end
out=fullfile(C.root,'analysis');if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,['va_step1_' name tag '.json']),'w');fprintf(fid,'%s',jsonencode(R,'PrettyPrint',true));fclose(fid);
end

function za=anchor(C,name,P) %#ok<INUSL>
M=make_model(name,P);za=M.spec.iface{C.anchor_iface}(P.x(P.in{2}));
end

% ---- identical to uq_rtm_20261117/scripts/score_stage2.m (local functions there) ----
function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=10;h=exp(-.5*((-30:30)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),30,1);y;repmat(y(end,:),30,1)],'valid');
y(P.z(P.in{1})<140,:)=0;
end

function s=score(a,r,z,x,W)
ix=x>=W{1}(1)&x<=W{1}(2);iz=find(z>=W{2}(1)&z<=W{2}(2));
rr=r(iz,ix);c0=cc(a(iz,ix),rr);cl=c0;
for L=[-2 -1 1 2],cl=max(cl,cc(a(iz+L,ix),rr));end
s=.5*(c0+cl);
end

function c=cc(a,b)
a=a(:)-mean(a(:));b=b(:)-mean(b(:));c=(a'*b)/sqrt((a'*a)*(b'*b));
end
