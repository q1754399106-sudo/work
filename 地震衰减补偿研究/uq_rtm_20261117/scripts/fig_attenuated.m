function fig_attenuated(name)
% Descriptive: reference image against the attenuated (uncompensated, -30 dB) image of model name at a
% common amplitude scale (no energy matching), image traces every 20 m, and wavenumber spectra of a
% shallow and a deep window. figures/attenuated_<name>.png
if nargin<1,name='C2';end
P=rtm_setup();S=load(fullfile(P.cache,['images2_' name '.mat']),'ref','m30');
z=P.z(P.in{1});x=P.x(P.in{2});ref=prep(S.ref.R0,P);att=prep(S.m30.R0,P);
iz=z>=300&z<=1550;xs=600:20:1800;[~,jx]=min(abs(x(:)-xs),[],1);sc=prctile(abs(ref(iz,jx)),99.5,'all');
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
f=figure('Visible','off','Position',[0 0 1800 760]);tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
pan={'参照（无衰减、无噪声）',ref;'衰减后（未补偿，道上 −30 dB），同一比例',att};
for p=1:2
 nexttile;hold on;A=pan{p,2};zz=z(iz)';
 for j=1:numel(jx)
  s=max(min(A(iz,jx(j))/sc*14,30),-30);x0=xs(j);
  fill([x0+max(s,0);x0*ones(size(s))],[zz;flipud(zz)],'k','EdgeColor','none');plot(x0+s,zz,'k','LineWidth',.3);
 end
 set(gca,'YDir','reverse');ylim([300 1550]);xlim([560 1840]);box on;grid on;title(pan{p,1});xlabel('x（m）');if p==1,ylabel('深度（m）');end
end
nexttile;hold on;cl=[0 0 0;.85 .4 .1];ws={[300 700],[1000 1500]};ls={'-','--'};
for w=1:2
 r=z>=ws{w}(1)&z<=ws{w}(2);n=1024;kz=(0:n/2)/(n*P.dx)*1000;              % cycles per km
 for p=1:2
  X=abs(fft(pan{p,2}(r,:).*hann(nnz(r)),n));s=mean(X(1:n/2+1,:),2);if p==1,s0=max(s);end
  plot(kz,20*log10(s/s0),'Color',cl(p,:),'LineStyle',ls{w},'LineWidth',1.4,'DisplayName',sprintf('%s，%d–%d m',strtok(pan{p,1},'（'),ws{w}(1),ws{w}(2)));
 end
end
xlim([0 40]);ylim([-45 3]);grid on;box on;xlabel('深度方向波数（周/km）');ylabel('振幅谱（dB，各窗相对参照最大值）');
legend('Location','southwest');title('振幅谱：浅窗（实线）与深窗（虚线）');
title(tl,sprintf('%s：参照与衰减后（未补偿）成像，同一振幅比例',name));
exportgraphics(f,fullfile(P.root,'figures',['attenuated_' name '.png']),'Resolution',140);close(f);
end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=10;h=exp(-.5*((-30:30)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),30,1);y;repmat(y(end,:),30,1)],'valid');
y(P.z(P.in{1})<140,:)=0;
end
