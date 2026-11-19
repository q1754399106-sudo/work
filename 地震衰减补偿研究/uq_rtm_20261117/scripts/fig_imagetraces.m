function fig_imagetraces(name)
% Descriptive: migrated-image traces (depth) of model name, data m30, for reference, uncompensated and
% compensated images (route A). Each image is scaled to the reference energy in 600-1550 m (display).
if nargin<1,name='C2';end
P=rtm_setup();S=load(fullfile(P.cache,['images2_' name '.mat']),'ref','m30','meta');rows=S.meta.combo.m30;
e1=find(strcmp(rows(:,1),'E1'));mm=[rows{e1,2}];z=P.z(P.in{1});x=P.x(P.in{2});iz=z>=600&z<=1550;
img=@(q,m) prep(S.m30.(q)(:,e1(mm==m)),P);ref=prep(S.ref.R0,P);
nrm=@(a) a*(norm(ref(iz,:),'fro')/max(norm(a(iz,:),'fro'),eps));
pan={'参照（无衰减、无噪声）',ref;'不补偿',nrm(prep(S.m30.R0,P));
 '反 Q 70 dB，m = 0.5',nrm(img('IQ70',.5));'W1（只用信噪比），m = 0.5',nrm(img('W1',.5));'UQ（本文），m = 0.5',nrm(img('UQW1',.5));
 '反 Q 70 dB，m = 1',nrm(img('IQ70',1));'W1，m = 1',nrm(img('W1',1));'UQ（本文），m = 1',nrm(img('UQW1',1))};
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
xs=600:20:1800;[~,jx]=min(abs(x(:)-xs),[],1);sc=prctile(abs(ref(iz,jx)),99.5,'all');
f=figure('Visible','off','Position',[0 0 1900 820]);tl=tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for p=1:8
 nexttile;hold on;A=pan{p,2};
 for j=1:numel(jx)
  s=max(min(A(iz,jx(j))/sc*14,30),-30);x0=xs(j);zz=z(iz)';
  fill([x0+max(s,0);x0*ones(size(s))],[zz;flipud(zz)],'k','EdgeColor','none');plot(x0+s,zz,'k','LineWidth',.3);
 end
 set(gca,'YDir','reverse');ylim([600 1550]);xlim([520 1880]);box on;grid on;title(pan{p,1});
 if mod(p,4)==1,ylabel('深度（m）');end;if p>4,xlabel('x（m）');end
end
title(tl,sprintf('%s 偏移成像的“道”（x = 600–1800 m，每 20 m 一道；道上 −30 dB；各图能量对齐到参照，只为显示）',name));
exportgraphics(f,fullfile(P.root,'figures',['imagetraces_' name '.png']),'Resolution',140);close(f);
% one trace at x = 1300 m, overlaid
[~,j]=min(abs(x-1300));zz=z(iz);
f=figure('Visible','off','Position',[0 0 1400 700]);tl=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;hold on
plot(zz,ref(iz,j),'k','LineWidth',1.4,'DisplayName','参照');plot(zz,pan{3,2}(iz,j),'Color',[.85 .55 .2],'DisplayName','反 Q 70 dB');
plot(zz,pan{4,2}(iz,j),'Color',[.35 .65 .35],'DisplayName','W1');plot(zz,pan{5,2}(iz,j),'Color',[.8 .15 .15],'LineWidth',1.3,'DisplayName','UQ（本文）');
grid on;box on;xlim([600 1550]);legend('Location','northeast');title('x = 1300 m，Q 偏小一半（m = 0.5，严重过补偿）');ylabel('成像振幅（对齐后）');
nexttile;hold on
plot(zz,ref(iz,j),'k','LineWidth',1.4,'DisplayName','参照');plot(zz,pan{2,2}(iz,j),'Color',[.6 .6 .6],'DisplayName','不补偿');
plot(zz,pan{7,2}(iz,j),'Color',[.35 .65 .35],'DisplayName','W1');plot(zz,pan{8,2}(iz,j),'Color',[.8 .15 .15],'LineWidth',1.3,'DisplayName','UQ（本文）');
grid on;box on;xlim([600 1550]);legend('Location','northeast');title('x = 1300 m，Q 准确（m = 1）');xlabel('深度（m）');ylabel('成像振幅（对齐后）');
exportgraphics(f,fullfile(P.root,'figures',['imagetrace1_' name '.png']),'Resolution',150);close(f);
end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=10;h=exp(-.5*((-30:30)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),30,1);y;repmat(y(end,:),30,1)],'valid');
y(P.z(P.in{1})<140,:)=0;
end
