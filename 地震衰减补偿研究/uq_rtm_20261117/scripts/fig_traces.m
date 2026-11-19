function fig_traces()
% Descriptive: wiggle traces of P1 shot 26 (x = 1200 m): elastic, viscoacoustic, viscoacoustic + noise
% (-30 dB), 3-55 Hz from the stored frequency samples; and the zero-offset trace overlaid.
P=rtm_setup();k=26;D=load(fullfile(P.cache,'data_P1.mat'),'ref','att','m30','C_m100');
t=P.t;E=exp(1i*2*pi*t*P.f);off=(P.ixr-P.ixs(k))*P.dx;
G={real(E*double(D.ref(:,:,k)).'),real(E*double(D.att(:,:,k)).'),real(E*double(D.m30(:,:,k)).'),real(E*double(D.C_m100(:,:,k)).')};
G{4}=G{4}*(norm(G{1}(:))/norm(G{4}(:)));                    % route C output has its own scale: match the reference energy (display)
nm={'无衰减（参照）','有衰减（Q 模型）','有衰减 + 道上 −30 dB 噪声','加噪声后逐道 UQ-W 补偿（m = 1）'};
sel=find(off>=0&off<=880&mod(off,80)==0);gn=(t/.5).^2;for q=1:4,G{q}=G{q}.*gn;end   % same t^2 display gain on all panels
sc=prctile(abs(G{1}(t>.6,sel)),99.5,'all');
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
f=figure('Visible','off','Position',[0 0 1800 760]);tl=tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
for q=1:4
 nexttile;hold on
 for j=1:numel(sel)
  s=max(min(G{q}(:,sel(j))/sc*35,70),-70);x0=off(sel(j));
  fill([x0+max(s,0);x0*ones(size(s))],[t;flipud(t)],'k','EdgeColor','none');
  plot(x0+s,t,'k','LineWidth',.5);
 end
 set(gca,'YDir','reverse');ylim([.25 1.6]);xlim([-60 940]);box on;grid on
 xlabel('偏移距（m）');if q==1,ylabel('时间（s）');end;title(nm{q});
end
title(tl,'P1 第 26 炮（x = 1200 m）的地震道（偏移距 0–880 m，每 80 m 一道；四幅乘同一 t² 增益、用同一比例，强振幅削顶）');
exportgraphics(f,fullfile(P.root,'figures','traces_P1.png'),'Resolution',150);close(f);
for q=1:4,G{q}=G{q}./gn;end
% zero-offset trace
j0=find(off==0);f=figure('Visible','off','Position',[0 0 1300 900]);tl=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;hold on
plot(t,G{1}(:,j0),'k','LineWidth',1.1,'DisplayName','无衰减');plot(t,G{2}(:,j0),'Color',[.85 .4 .1],'LineWidth',1.1,'DisplayName','有衰减');
xlim([.25 1.6]);grid on;box on;legend('Location','northeast');title('零偏移距道：无衰减与有衰减（同一比例）');ylabel('压力（相对）');
nexttile;hold on
plot(t,G{3}(:,j0),'Color',[.5 .5 .5],'LineWidth',.8,'DisplayName','有衰减 + −30 dB 噪声');plot(t,G{2}(:,j0),'Color',[.85 .4 .1],'LineWidth',1.1,'DisplayName','有衰减（无噪声）');
xlim([.25 1.6]);grid on;box on;legend('Location','northeast');xlabel('时间（s）');ylabel('压力（相对）');
title('同一道加噪声后：0.8 s 以后的反射已低于噪声');
nexttile;hold on
plot(t,G{1}(:,j0),'k','LineWidth',1.1,'DisplayName','无衰减（目标）');plot(t,G{4}(:,j0),'Color',[.2 .45 .8],'LineWidth',1.1,'DisplayName','加噪声后逐道 UQ-W 补偿（m = 1，能量对齐）');
xlim([.25 1.6]);grid on;box on;legend('Location','northeast');xlabel('时间（s）');ylabel('压力（相对）');
title('逐道补偿：浅层波形变尖、接近目标；深部噪声主导，补不回来');
exportgraphics(f,fullfile(P.root,'figures','trace0_P1.png'),'Resolution',150);close(f);
end
