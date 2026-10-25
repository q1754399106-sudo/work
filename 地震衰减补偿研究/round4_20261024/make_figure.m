function make_figure()
% Oracle vs data-only layered/global ratio per scene.
root=fileparts(mfilename('fullpath'));o=[jsondecode(fileread(fullfile(root,'analysis','oracle_V4_F1_F2.json')));jsondecode(fileread(fullfile(root,'analysis','oracle_F3_F4.json')))];
d=[jsondecode(fileread(fullfile(root,'analysis','datalayer_V4_F1_F2.json')));jsondecode(fileread(fullfile(root,'analysis','datalayer_F3_F4.json')))];
lab=arrayfun(@(r)sprintf('%s Q%+d',r.set,r.Q_percent),o,'UniformOutput',false);
f=figure('Visible','off','Position',[0 0 1400 520]);b=bar([[o.layered_vs_global]',[d.ratio]']);
b(1).FaceColor=[.25 .55 .85];b(2).FaceColor=[.9 .45 .2];yline(1,'k-');yline(.97,'k:');
set(gca,'XTick',1:numel(lab),'XTickLabel',lab,'XTickLabelRotation',40);ylabel('layered error / global error');ylim([.7 1.4]);grid on;
legend({'oracle (truth-chosen layer factors)','data only (lasso objective)','no change','decision line 0.97'},'Location','northwest');
title(sprintf('Depth-segmented QScan: oracle mean %.3f vs data-only mean %.3f',mean([o.layered_vs_global]),mean([d.ratio])));
exportgraphics(f,fullfile(root,'analysis','layered_vs_global.png'),'Resolution',160);close(f);
end
