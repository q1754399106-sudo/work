function fld_phaserot(outdir)
% Second post-hoc check of stage 4 (PLAN_PHASEROT_ZH.md): constant phase rotation of the F-5 well-side trace.
% Tie exactly as fld_eval.m (40 m detrend, +-2 traces, R0 -10 dB depth-wavenumber band, +-3 sample shift), plus a
% constant rotation tr_th = cos(th) tr - sin(th) H(tr), th = -180:2:178 deg. Image sets: original (f_ref 100), A
% (f_ref 25), B (amplitude only). Scale x_bar = median two-way t* (model Q 120) over F-5 reservoir points and shots.
% Run from the main checkout's scripts folder (cache there); writes <outdir>/stage4_phaserot.json.
P=fld_setup();z=P.z(P.in{1});x=P.x(P.in{2});
sets={'orig','images.mat';'A','images_fref25.mat';'B','images_amp.mat'};meth={'PH','UQW1','IQ70','W1'};
th=(-180:2:178)*pi/180;f0=14;
% R0 band, as fld_eval
I=load(fullfile(P.cache,'images.mat'));ms=I.meta.m;
wz=z>=2600&z<=3200;wx=x>=-2500&x<=1500;n=256;kz=(0:n/2)/(n*P.dx)*1000;
A0=prep(I.R0(:,1),P);s0=mean(abs(fft(A0(wz,wx).*hann(nnz(wz)),n)),2);s0=s0(1:n/2+1);b=band(kz,s0);lo=b(1);hi=b(2);
% F-5 well
fid=fopen(fullfile(P.cache,'well_F-5.bin'));a=fread(fid,Inf,'single');fclose(fid);a=reshape(a,4,[])';
zw=a(:,1);sw=a(:,2);rw=a(:,4);[~,iz]=min(abs(z(:)-zw'),[],1);rb=bandpass_z(rw,P.dx,lo,hi);
jx=zeros(numel(zw),1);for i=1:numel(zw),[~,jx(i)]=min(abs(x-sw(i)));end
% t* scale
T=load(fullfile(P.cache,'ts.mat'),'TS');lin=sub2ind([P.nzi P.nxi],iz(:),jx(:));xb=median(double(T.TS(lin,:)),'all');clear T
R=struct('f0',f0,'xbar',xb,'theta0_deg_m1',2*f0*xb*180/pi,'m',ms,'band',[lo hi],'theta_deg',th*180/pi);
[R.R0.c0,R.R0.cstar,R.R0.thstar]=tie(A0,iz,jx,rb,P,lo,hi,th);
for s=1:size(sets,1)
 if s>1,I=load(fullfile(P.cache,sets{s,2}));end
 for q=meth
  for k=1:numel(ms)
   [c0,cs,ts]=tie(prep(I.(q{1})(:,k),P),iz,jx,rb,P,lo,hi,th);
   R.(sets{s,1}).(q{1}).c0(k)=c0;R.(sets{s,1}).(q{1}).cstar(k)=cs;R.(sets{s,1}).(q{1}).thstar(k)=ts;
  end
 end
end
fid=fopen(fullfile(outdir,'stage4_phaserot.json'),'w');fprintf(fid,'%s',jsonencode(R,'PrettyPrint',true));fclose(fid);
fprintf('x_bar %.4f s, theta0(m=1) %.1f deg, R0 c0 %.3f c* %.3f th* %d\n',xb,R.theta0_deg_m1,R.R0.c0,R.R0.cstar,R.R0.thstar);
end

function [c0,cs,ts]=tie(A,iz,jx,rb,P,lo,hi,th)
tr=zeros(numel(iz),1);for i=1:numel(iz),j=jx(i);tr(i)=mean(A(iz(i),max(j-2,1):min(j+2,size(A,2))));end
tr=bandpass_z(tr,P.dx,lo,hi);h=hilb(tr);c0=best(tr,rb);cs=-Inf;ts=0;
for k=1:numel(th),c=best(cos(th(k))*tr-sin(th(k))*h,rb);if c>cs,cs=c;ts=round(th(k)*180/pi);end,end
end

function c=best(tr,rb),c=-Inf;for L=-3:3,c=max(c,cc(circshift(tr,L),rb));end,end

function h=hilb(y)
n=2^nextpow2(4*numel(y));Y=fft(y,n);w=zeros(n,1);w(2:n/2)=-1i;w(n/2+2:n)=1i;h=real(ifft(Y.*w));h=h(1:numel(y));
end

function b=band(k,s)
s=s(:)';[m,i]=max(s);t=m*10^(-10/20);l=find(s(1:i)<t,1,'last');h=i-1+find(s(i:end)<t,1);
if isempty(l),l=1;end;if isempty(h),h=numel(s);end;b=[k(l) k(h)];
end

function y=bandpass_z(r,dx,lo,hi)
n=2^nextpow2(4*numel(r));R=fft(r,n);k=abs([0:n/2 -n/2+1:-1]')/(n*dx)*1000;y=real(ifft(R.*(k>=lo&k<=hi)));y=y(1:numel(r));
end

function c=cc(a,b),a=a-mean(a);b=b-mean(b);c=(a'*b)/sqrt((a'*a)*(b'*b));end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=4;h=exp(-.5*((-12:12)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),12,1);y;repmat(y(end,:),12,1)],'valid');y(P.z(P.in{1})<150,:)=0;
end
