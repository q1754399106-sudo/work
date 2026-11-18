function fig_stage2()
% Descriptive figures of stage 2 (development, D1-D4): S-m curves, UQ minus best feasible baseline,
% image panels of D2 (m30) at m = 0.5 and m = 1. Output figures/stage2_*.png.
P=rtm_setup();out=fullfile(P.root,'figures');if ~isfolder(out),mkdir(out);end
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
names={'D1','D2','D3','D4'};desc={'D1 断层+两透镜','D2 背斜+气烟囱','D3 Q 横向渐变','D4 浅层低 Q 薄层'};
meth={'R0','PH','IQ70','LP','W0','UQW0'};lab={'不补偿','只校正相位','反 Q（70 dB）','完全补偿+低通','W0（只用信噪比）','UQ（本文）'};
col=[.6 .6 .6;.30 .55 .80;.85 .55 .20;.55 .40 .70;.35 .65 .35;.80 .15 .15];
R=cell(1,4);for k=1:4,R{k}=jsondecode(fileread(fullfile(P.root,'analysis',['stage2_' names{k} '.json'])));end
% 1) S-m curves, m30 and m40, deep window
for ds={'m30','m40'}
 f=figure('Visible','off','Position',[0 0 1300 400]);tl=tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
 for k=1:4
  T=R{k}.(ds{1}).deep.E1;ax=nexttile;hold on
  for q=1:numel(meth)
   y=T.(meth{q});lw=1.2;if strcmp(meth{q},'UQW0'),lw=2.6;end
   plot(T.m,y,'-o','Color',col(q,:),'LineWidth',lw,'MarkerSize',3,'MarkerFaceColor',col(q,:),'DisplayName',lab{q});
  end
  xline(1,':k','HandleVisibility','off');set(gca,'XScale','log','XTick',T.m,'XTickLabel',compose('%g',T.m));
  xlim([.48 2.1]);ylim([0 1]);grid on;box on;title(desc{k});xlabel('Q 乘子 m（<1：所用 Q 偏小，过补偿）');
  if k==1,ylabel('S（与参照成像的相似度）');end
  if k==4,lg=legend('NumColumns',6,'FontSize',9);lg.Layout.Tile='south';end
 end
 title(tl,sprintf('各方法随 Q 误差的成像相似度（深部窗，道上噪声 %s dB）',strrep(strrep(ds{1},'m3','-3'),'m4','-4')));
 exportgraphics(f,fullfile(out,['stage2_S_m_' ds{1} '.png']),'Resolution',160);close(f);
end
% 2) UQ minus best feasible baseline
f=figure('Visible','off','Position',[0 0 1100 380]);mk={'o','s','^','d'};cl=lines(4);
for a=1:2
 ds={'m30','m40'};subplot(1,2,a);hold on
 for k=1:4
  T=R{k}.(ds{a}).deep.E1;plot(T.m,T.d_Bstar,['-' mk{k}],'Color',cl(k,:),'LineWidth',1.6,'MarkerFaceColor',cl(k,:),'DisplayName',desc{k});
 end
 yline(0,'k','HandleVisibility','off');yline(-.05,'--','Color',[.5 .5 .5],'HandleVisibility','off');yline(.05,'--','Color',[.5 .5 .5],'HandleVisibility','off');
 xline(1,':k','HandleVisibility','off');set(gca,'XScale','log','XTick',[.5 .6 .7 .85 1 1.2 1.4 1.7 2],'XTickLabel',compose('%g',[.5 .6 .7 .85 1 1.2 1.4 1.7 2]));
 xlim([.48 2.1]);ylim([-.15 .3]);grid on;box on;xlabel('Q 乘子 m');ylabel('UQ − 最好可行基线（S）');
 title(sprintf('道上噪声 %s dB，深部窗',strrep(strrep(ds{a},'m3','-3'),'m4','-4')));if a==1,legend('Location','northeast');end
end
sgtitle('UQ 相对“逐情景最好的可行基线”的差值（虚线：草案阈值 ±0.05）');
exportgraphics(f,fullfile(out,'stage2_delta.png'),'Resolution',160);close(f);
% 3) image panels: D2, m30, deep part, m = 0.5 and 1
S=load(fullfile(P.cache,'images2_D2.mat'),'ref','m30','meta');rows=S.meta.combo.m30;
e1=find(strcmp(rows(:,1),'E1'));mm=[rows{e1,2}];z=P.z(P.in{1});x=P.x(P.in{2});iz=z>=600;
pick=@(q,m) S.m30.(q)(:,e1(mm==m));ref=prep(S.ref.R0,P);c=prctile(abs(ref(z>=820,:)),99,'all');
pan={'参照（无衰减、无噪声）',ref;'不补偿',prep(S.m30.R0,P)};
for m=[.5 1]
 for q={'IQ70','W0','UQW0'}
  pan(end+1,:)={sprintf('%s，m = %g',lab{strcmp(meth,q{1})},m),prep(pick(q{1},m),P)}; %#ok<AGROW>
 end
end
f=figure('Visible','off','Position',[0 0 1500 470]);ord=[1 3 4 5 2 6 7 8];tl=tiledlayout(2,4,'TileSpacing','tight','Padding','compact');
for k=1:8
 nexttile;img=pan{ord(k),2};
 if ord(k)>1,img=img*(norm(ref(iz,:),'fro')/max(norm(img(iz,:),'fro'),eps));end   % common scale for display only
 imagesc(x,z(iz),img(iz,:),[-c c]);colormap(gray);axis image;title(pan{ord(k),1},'FontSize',9);
 if mod(k,4)==1,ylabel('深度（m）');end
 if k>4,xlabel('x（m）');end
 rectangle('Position',[1050 850 300 450],'EdgeColor',[.9 .2 .2],'LineWidth',1);
end
title(tl,'D2 成像对比（道上 −30 dB；上排 m = 0.5 严重过补偿，下排 m = 1 Q 准确；红框为气烟囱下方窗）');
exportgraphics(f,fullfile(out,'stage2_images_D2.png'),'Resolution',160);close(f);
end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=10;h=exp(-.5*((-30:30)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),30,1);y;repmat(y(end,:),30,1)],'valid');
y(P.z(P.in{1})<140,:)=0;
end
