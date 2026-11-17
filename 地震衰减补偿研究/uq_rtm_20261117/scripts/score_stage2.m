function R=score_stage2(name,vscale)
% Scores of PLAN_STAGE2_ZH.md section 6 for one model (and velocity variant). Writes
% analysis/stage2_<name><tag>.json.
if nargin<2,vscale=1;end
P=rtm_setup();C=s2_cfg(name,vscale);I=load(fullfile(P.cache,['images2_' name C.tag '.mat']));
z=P.z(P.in{1});x=P.x(P.in{2});ref=prep(I.ref.R0,P);wn=fieldnames(C.win);
base={'R0','PH','IQ30','IQ70','LP','W0'};meth=[{'R0'} C.methods];
R=struct('model',name,'tag',C.tag,'fc',I.meta.fc,'win',C.win);
for d=1:size(C.sets,1)
 ds=C.sets{d,1};rows=I.meta.combo.(ds);
 for w=1:numel(wn)
  W=C.win.(wn{w});S=nan(numel(meth),size(rows,1));
  for q=1:numel(meth)
   for r=1:size(rows,1)
    if strcmp(meth{q},'R0'),img=I.(ds).R0;else,img=I.(ds).(meth{q})(:,r);end
    S(q,r)=score(prep(img,P),ref,z,x,W);
   end
  end
  forms=unique(rows(:,1),'stable');
  for fo=forms'
   k=strcmp(rows(:,1),fo{1});T=struct('m',[rows{k,2}]);
   for q=1:numel(meth),T.(meth{q})=S(q,k);end
   [~,ib]=ismember(base,meth);T.Bstar=max(S(ib,k),[],1);u=S(strcmp(meth,'UQW0'),k);
   T.d_Bstar=u-T.Bstar;T.d_W0=u-S(strcmp(meth,'W0'),k);T.d_LP=u-S(strcmp(meth,'LP'),k);
   ms=T.m;
   if any(ms<1),T.M1=mean(T.d_Bstar(ms<1));end
   if any(ms==1),T.M2=T.d_Bstar(ms==1);end
   if any(ms>1),T.M3=min(T.d_Bstar(ms>1));end
   T.M4=struct();for b=base,T.M4.(b{1})=min(u)-min(S(strcmp(meth,b{1}),k));end
   R.(ds).(wn{w}).(fo{1})=T;
  end
 end
end
out=fullfile(P.root,'analysis');if ~isfolder(out),mkdir(out);end
fid=fopen(fullfile(out,['stage2_' name C.tag '.json']),'w');fprintf(fid,'%s',jsonencode(R,'PrettyPrint',true));fclose(fid);
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
