function fld_eval(imgfile,sfx)
if nargin<1,imgfile='images.mat';end
if nargin<2,sfx='';end
% Stage-4 descriptive indicators (PLAN_STAGE4_ZH.md section 4) and figures. analysis/stage4_field.json,
% figures/field_*.png.
P=fld_setup();I=load(fullfile(P.cache,imgfile));S=load(fullfile(P.cache,'inv.mat'));
z=P.z(P.in{1});x=P.x(P.in{2});meth=[{'R0'} I.meta.methods];ms=I.meta.m;
img=@(q,k) prep(I.(q)(:,min(k,size(I.(q),2))),P);
wz=z>=2600&z<=3200;wx=x>=-2500&x<=1500;
% ---- 1. bandwidth in the reservoir window ----
n=256;kz=(0:n/2)/(n*P.dx)*1000;                                             % cycles / km
spec=@(A) mean(abs(fft(A(wz,wx).*hann(nnz(wz)),n)),2);
s0=spec(img('R0',1));s0=s0(1:n/2+1);R.kz=kz;
bw=@(s) band(kz,s);[R.R0.bw10,R.R0.cen]=bw(s0);lo=R.R0.bw10(1);hi=R.R0.bw10(2);
for q=meth(2:end)
 for k=1:numel(ms)
  s=spec(img(q{1},k));s=s(1:n/2+1);[b,c]=bw(s);R.(q{1}).bw10(k,:)=b;R.(q{1}).cen(k)=c;
 end
end
% ---- 2. well ties ----
for w={'F-5','F-14'}
 fid=fopen(fullfile(P.cache,['well_' w{1} '.bin']));a=fread(fid,Inf,'single');fclose(fid);a=reshape(a,4,[])';
 zw=a(:,1);sw=a(:,2);rw=a(:,4);[~,iz]=min(abs(z(:)-zw'),[],1);
 rb=bandpass_z(rw,P.dx,lo,hi);                                               % R0 -10 dB band
 for q=meth
  for k=1:(1+(numel(ms)-1)*~strcmp(q{1},'R0'))
   A=img(q{1},k);tr=zeros(numel(zw),1);
   for i=1:numel(zw),[~,j]=min(abs(x-sw(i)));tr(i)=mean(A(iz(i),max(j-2,1):min(j+2,numel(x))));end
   tr=bandpass_z(tr,P.dx,lo,hi);best=-Inf;for L=-3:3,c=cc(circshift(tr,L),rb);if c>best,best=c;end,end   % image trace in the same band (plan 4.2; fixed 11-19)
   R.(q{1}).(['tie_' strrep(w{1},'-','')])(k)=best;
  end
 end
 R.(['well_' strrep(w{1},'-','')])=struct('z',[min(zw) max(zw)],'s',[min(sw) max(sw)]);
end
% ---- 3. image-domain SNR in the reservoir window ----
wr=wz'&wx;R.snr_med=1./median(S.iv(wr(:),:),1);R.f=P.f;R.fcut=I.meta.fcut;R.m=ms;R.methods=meth;
out=fullfile(P.root,'analysis');fid=fopen(fullfile(out,['stage4_field' sfx '.json']),'w');fprintf(fid,'%s',jsonencode(R,'PrettyPrint',true));fclose(fid);
% ---- figures ----
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
zz=z>=1500&z<=3400;ref=img('R0',1);
pans={'R0',1,'不补偿';'IQ70',2,'反 Q 70 dB，m = 0.6';'LP',2,'完全补偿+低通，m = 0.6';'W1',2,'W1，m = 0.6';'UQW1',2,'UQ，m = 0.6';
      'IQ70',5,'反 Q 70 dB，m = 1';'W1',5,'W1，m = 1';'UQW1',5,'UQ，m = 1';'UQW1',7,'UQ，m = 1.4'};
f=figure('Visible','off','Position',[0 0 1800 1000]);tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for p=1:size(pans,1)
 A=img(pans{p,1},pans{p,2});A=A*(norm(ref(zz,:),'fro')/max(norm(A(zz,:),'fro'),eps));c=prctile(abs(ref(zz,:)),99,'all');
 nexttile;imagesc(x/1000,z(zz),A(zz,:),[-c c]);colormap(gray);title(pans{p,3});hold on
 plot([-1.275 -1.118],[2570 3170],'r-','LineWidth',1.5);if mod(p,3)==1,ylabel('深度（m）');end;if p>6,xlabel('测线坐标（km，平台为 0）');end
end
exportgraphics(f,fullfile(P.root,'figures',['field_images' sfx '.png']),'Resolution',120);close(f);
f=figure('Visible','off','Position',[0 0 1500 450]);tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
cl=lines(8);mm={'PH','IQ30','IQ70','LP','W0','W1','UQW0','UQW1'};
nexttile;hold on;for i=1:numel(mm),plot(ms,R.(mm{i}).cen,'-o','Color',cl(i,:),'LineWidth',1+1.5*strcmp(mm{i},'UQW1'),'DisplayName',mm{i});end
yline(R.R0.cen,'k--','R0','HandleVisibility','off');set(gca,'XScale','log');grid on;xlabel('Q 乘子 m（名义 Q = 120）');ylabel('谱质心（周/km）');title('储层窗深度波数谱质心');legend('Location','best','NumColumns',2)
nexttile;hold on;for i=1:numel(mm),plot(ms,R.(mm{i}).tie_F5,'-o','Color',cl(i,:),'LineWidth',1+1.5*strcmp(mm{i},'UQW1'));end
yline(R.R0.tie_F5,'k--','R0');set(gca,'XScale','log');grid on;xlabel('Q 乘子 m');ylabel('相关系数');title('F-5 井旁相关');
nexttile;plot(P.f,10*log10(R.snr_med),'k','LineWidth',1.4);grid on;xlabel('频率（Hz）');ylabel('dB');title('储层窗成像域信噪比（估计，中位数）');
exportgraphics(f,fullfile(P.root,'figures',['field_metrics' sfx '.png']),'Resolution',140);close(f);
end

function [b,c]=band(k,s)
s=s(:)';[m,i]=max(s);t=m*10^(-10/20);l=find(s(1:i)<t,1,'last');h=i-1+find(s(i:end)<t,1);
if isempty(l),l=1;end;if isempty(h),h=numel(s);end;b=[k(l) k(h)];c=sum(k.*s.^2)/sum(s.^2);
end

function y=bandpass_z(r,dx,lo,hi)
n=2^nextpow2(4*numel(r));R=fft(r,n);k=abs([0:n/2 -n/2+1:-1]')/(n*dx)*1000;y=real(ifft(R.*(k>=lo&k<=hi)));y=y(1:numel(r));
end

function c=cc(a,b),a=a-mean(a);b=b-mean(b);c=(a'*b)/sqrt((a'*a)*(b'*b));end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=4;h=exp(-.5*((-12:12)/s).^2);h=h/sum(h);   % 40 m smoothing at 10 m grid
y=y-conv2(h',1,[repmat(y(1,:),12,1);y;repmat(y(end,:),12,1)],'valid');y(P.z(P.in{1})<150,:)=0;
end
