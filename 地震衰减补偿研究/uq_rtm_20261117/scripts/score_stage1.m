function score_stage1(name)
% Scores of PLAN_STAGE1_ZH.md (run only after the plan is committed). Writes analysis/stage1_<name>.json.
P=rtm_setup();I=load(fullfile(P.cache,['images_' name '.mat']));
z=P.z(P.in{1});x=P.x(P.in{2});
win=struct('sh',{{[300 2100],[150 400]}},'mid',{{[300 2100],[420 800]}},'deep',{{[300 2100],[820 1550]}},'lens',{{[600 1200],[700 1100]}});
ref=prep(I.ref.R0,P);wn=fieldnames(win);
R=struct('model',name,'windows',win);
meth=[{'R0'} I.meta.methods];base={'R0','PH','IQ30','IQ70','W0'};iqs={'IQ20','IQ30','IQ40','IQ50','IQ70'};
for d=1:size(I.meta.m_full,1)
 ds=I.meta.m_full{d,1};ms=I.meta.m_full{d,2};
 for w=1:numel(wn)
  W=win.(wn{w});S=nan(numel(meth),numel(ms));
  for q=1:numel(meth)
   for k=1:numel(ms)
    if strcmp(meth{q},'R0'),img=I.(ds).R0;else,img=I.(ds).(meth{q})(:,k);end
    S(q,k)=score(prep(img,P),ref,z,x,W);
   end
  end
  T=struct('m',ms);for q=1:numel(meth),T.(meth{q})=S(q,:);end
  [~,ib]=ismember(base,meth);[~,ii]=ismember(iqs,meth);
  T.Bstar=max(S(ib,:),[],1);T.IQoracle=max(S(ii,:),[],1);
  u=S(strcmp(meth,'UQW0'),:);T.d_Bstar=u-T.Bstar;T.d_W0=u-S(strcmp(meth,'W0'),:);T.d_IQoracle=u-T.IQoracle;
  lo=ms<1;hi=ms>1;one=ms==1;
  T.M1_draft=mean(T.d_Bstar(lo));T.M2_draft=T.d_Bstar(one);T.M3_draft=T.d_Bstar(hi);
  T.M4_draft=struct();for b=base,T.M4_draft.(b{1})=min(u)-min(S(strcmp(meth,b{1}),:));end
  R.(ds).(wn{w})=T;
 end
 if isfield(I.(ds),'inv_med'),R.(ds).inv_med=I.(ds).inv_med;end
end
% route C against route A (UQW0) on the route-C dataset
dsC=P.routeC.data;msA=I.meta.m_full{strcmp(I.meta.m_full(:,1),dsC),2};
for w=1:numel(wn)
 W=win.(wn{w});C=struct('m',P.routeC.m);
 for k=1:numel(P.routeC.m)
  m=P.routeC.m(k);tag=sprintf('C_m%03d',round(100*m));
  C.S_C(k)=score(prep(I.(tag).R0,P),ref,z,x,W);C.S_A(k)=score(prep(I.(dsC).UQW0(:,msA==m),P),ref,z,x,W);
 end
 C.d_AC=C.S_A-C.S_C;C.mean_d_AC=mean(C.d_AC);R.routeC.(wn{w})=C;
end
R.route_choice=ternary(R.routeC.deep.mean_d_AC>=-.01,'A','C');
out=fullfile(P.root,'analysis');if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,['stage1_' name '.json']),'w');fprintf(fid,'%s',jsonencode(R,'PrettyPrint',true));fclose(fid);
fprintf('route choice: %s (mean d_AC deep %.4f)\n',R.route_choice,R.routeC.deep.mean_d_AC);
end

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

function o=ternary(c,a,b)
if c,o=a;else,o=b;end
end
