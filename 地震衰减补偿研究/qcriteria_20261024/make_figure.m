function make_figure()
root=fileparts(mfilename('fullpath'));r=jsondecode(fileread(fullfile(root,'analysis','criteria.json')));
crit={'C0_QScan','C1_evidence','C2_GCV','C3_holdout','C4_spectral'};data={'M1_full','M5_blocky_sparse','M1_noisy'};B=zeros(5,3);
for i=1:5,for k=1:3,s=strcmp({r.criterion},crit{i})&strcmp({r.data},data{k});B(i,k)=mean(abs([r(s).bias_log]));end,end
f=figure('Visible','off','Position',[0 0 1100 480]);h=bar(B);h(1).FaceColor=[.2 .4 .8];h(2).FaceColor=[.85 .3 .3];h(3).FaceColor=[.5 .5 .5];
yline(.054,'k--','success line (1 grid step)');set(gca,'XTickLabel',strrep(crit,'_',' '));ylabel('mean |ln(c_{sel}/c_{opt})| over F1-F4');grid on;
legend({'M1 dense, noise-free','M5 sparse, noise-free','M1 dense, noisy (archived)'},'Location','northwest');
title('Q-selection criteria: none reaches the registered success line');
exportgraphics(f,fullfile(root,'analysis','criteria_bias.png'),'Resolution',160);close(f);
fprintf('%-12s %6s %6s %6s\n','criterion',data{:});for i=1:5,fprintf('%-12s %6.3f %6.3f %6.3f\n',crit{i},B(i,:));end
end
