function make_figure()
root=fileparts(mfilename('fullpath'));a=[jsondecode(fileread(fullfile(root,'analysis','variants_F1_F2.json')));jsondecode(fileread(fullfile(root,'analysis','variants_F3_F4.json')))];
b=jsondecode(fileread(fullfile(root,'analysis','variants_F1_F2_F3_F4_M5_blocky_sparse.json')));r=[a;b];
sets={'F1','F2','F3','F4'};vars={'M1_full','M2_primaries','M3_noloss','M4_noloss_macroQ','M5_blocky_sparse'};B=zeros(4,5);
for i=1:4,for k=1:5,B(i,k)=r(strcmp({r.set},sets{i})&strcmp({r.variant},vars{k})).bias_log;end,end
f=figure('Visible','off','Position',[0 0 1200 520]);h=bar(B);hold on;
cols=[.2 .4 .8;.35 .65 .9;.5 .8 .6;.95 .7 .3;.85 .3 .3];for k=1:5,h(k).FaceColor=cols(k,:);end
yline(log(1.28),'k:','ln 1.28');set(gca,'XTickLabel',sets);ylabel('selection bias  ln(c_{sel} / c_{opt})');grid on;
legend({'M1 full (archived)','M2 no internal multiples','M3 + no transmission loss','M4 + smooth Q','M5 + blocky sparse reflectivity (post hoc)'},'Location','northwest');
title(sprintf('QScan bias on noise-free Kennett data: mean M1 %.3f, M4 %.3f, M5 %.3f',mean(B(:,1)),mean(B(:,4)),mean(B(:,5))));
exportgraphics(f,fullfile(root,'analysis','bias_by_variant.png'),'Resolution',160);close(f);
end
