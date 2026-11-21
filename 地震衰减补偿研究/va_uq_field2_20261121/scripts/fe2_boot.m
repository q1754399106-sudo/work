function B=fe2_boot(tag)
% Bootstrap intervals for the F-5 tie of fld_eval.m (PLAN_ZH.md section 6). fld_eval.m has no interval, so
% this script recomputes its F-5 tie with the same code (prep, trace extraction +-2 traces, R0 -10 dB box
% band, bandpass_z, cc, shift -3..3 samples, copied verbatim) and checks every point value against
% analysis/fe2_eval<tag>.json (|diff| <= 1e-9). Then, at each method's own best shift (as
% field_application_20261028), a circular moving-block bootstrap over the depth samples of the tie
% window: block 10 samples (100 m), 2000 resamples, seed 20261121, same resampled indices for all
% methods; 90% percentile intervals of the correlation and of (method - R0). Descriptive only.
if nargin<1,tag='';end
[C,P]=fe2_cfg();I=load(fullfile(C.cache,['images_fe2' tag '.mat']));
R=jsondecode(fileread(fullfile(C.root,'analysis',['fe2_eval' tag '.json'])));
z=P.z(P.in{1});x=P.x(P.in{2});lo=R.R0.bw10(1);hi=R.R0.bw10(2);ms=I.meta.m;meth=[{'R0'} I.meta.methods(:)'];
fid=fopen(fullfile(C.s4cache,'well_F-5.bin'));a=fread(fid,Inf,'single');fclose(fid);a=reshape(a,4,[])';
zw=a(:,1);sw=a(:,2);rw=a(:,4);[~,iz]=min(abs(z(:)-zw'),[],1);rb=bandpass_z(rw,P.dx,lo,hi);n=numel(rb);
TR=struct();maxdiff=0;
for q=meth
 K=1+(numel(ms)-1)*~strcmp(q{1},'R0');TR.(q{1})=zeros(n,K);
 for k=1:K
  A=prep(I.(q{1})(:,min(k,size(I.(q{1}),2))),P);tr=zeros(n,1);
  for i=1:n,[~,j]=min(abs(x-sw(i)));tr(i)=mean(A(iz(i),max(j-2,1):min(j+2,numel(x))));end
  tr=bandpass_z(tr,P.dx,lo,hi);best=-Inf;bl=0;for L=-3:3,c=cc(circshift(tr,L),rb);if c>best,best=c;bl=L;end,end
  TR.(q{1})(:,k)=circshift(tr,bl);maxdiff=max(maxdiff,abs(best-R.(q{1}).tie_F5(k)));
 end
end
fprintf('point values vs fld_eval json: max |diff| = %.3g\n',maxdiff);if maxdiff>1e-9,error('tie recomputation differs from fld_eval');end
nbk=10;nrep=2000;rng(20261121,'twister');nb=ceil(n/nbk);IX=zeros(n,nrep);
for r=1:nrep,st=randi(n,nb,1);ix=mod(st+(0:nbk-1)-1,n)+1;ix=reshape(ix',[],1);IX(:,r)=ix(1:n);end
c0=zeros(nrep,1);for r=1:nrep,c0(r)=cc(TR.R0(IX(:,r)),rb(IX(:,r)));end
B=struct('block',nbk,'nrep',nrep,'seed',20261121,'n',n,'m',ms,'point_check_maxdiff',maxdiff);
B.R0=struct('tie',R.R0.tie_F5,'ci90',prctile(c0,[5 95]));
for q=meth(2:end)
 K=size(TR.(q{1}),2);ci=zeros(K,2);cd=ci;
 for k=1:K
  cb=zeros(nrep,1);for r=1:nrep,cb(r)=cc(TR.(q{1})(IX(:,r),k),rb(IX(:,r)));end
  ci(k,:)=prctile(cb,[5 95]);cd(k,:)=prctile(cb-c0,[5 95]);
 end
 B.(q{1})=struct('tie',R.(q{1}).tie_F5(:)','ci90',ci,'diff_R0',R.(q{1}).tie_F5(:)'-R.R0.tie_F5,'diff_ci90',cd);
end
fid=fopen(fullfile(C.root,'analysis',['fe2_boot' tag '.json']),'w');fprintf(fid,'%s',jsonencode(B,'PrettyPrint',true));fclose(fid);
fprintf('fe2_boot%s written %s\n',tag,datestr(now));
end

function y=bandpass_z(r,dx,lo,hi)
n=2^nextpow2(4*numel(r));R=fft(r,n);k=abs([0:n/2 -n/2+1:-1]')/(n*dx)*1000;y=real(ifft(R.*(k>=lo&k<=hi)));y=y(1:numel(r));
end

function c=cc(a,b),a=a-mean(a);b=b-mean(b);c=(a'*b)/sqrt((a'*a)*(b'*b));end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=4;h=exp(-.5*((-12:12)/s).^2);h=h/sum(h);   % 40 m smoothing at 10 m grid
y=y-conv2(h',1,[repmat(y(1,:),12,1);y;repmat(y(end,:),12,1)],'valid');y(P.z(P.in{1})<150,:)=0;
end
