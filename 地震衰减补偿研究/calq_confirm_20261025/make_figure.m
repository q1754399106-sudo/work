function make_figure()
root=fileparts(mfilename('fullpath'));M=readtable(fullfile(root,'analysis','confirm_metrics.csv'));Q=jsondecode(fileread(fullfile(root,'analysis','confirm_q.json')));
f=figure('Visible','off','Position',[0 0 1500 500]);tiledlayout(1,2,'TileSpacing','compact');
nexttile;keys={};V=[];
for m={'C1','C2','C3'},for c={'Moderate','Strong'}
 s=strcmp(M.model,m{1})&strcmp(M.condition,c{1})&strcmp(M.method,'Standard');
 V(end+1,:)=[M.overall(s&strcmp(M.calibration,'current')),M.overall(s&strcmp(M.calibration,'SharedGhost'))];keys{end+1}=sprintf('%s %s',m{1},c{1}); %#ok<AGROW>
end,end
b=bar(V);b(1).FaceColor=[.2 .45 .8];b(2).FaceColor=[.9 .6 .1];set(gca,'XTickLabel',keys,'XTickLabelRotation',25);ylabel('Standard overall relative error');grid on;ylim([0 .8]);
legend({'current calibration','SharedGhost'},'Location','northwest');title('Exploratory confirmation: section error');
nexttile;B=[];lab={};
for m={'C1','C2','C3'},for c={'current','SharedGhost'}
 q=Q(strcmp({Q.model},m{1})&strcmp({Q.calibration},c{1}));B(end+1,:)=q.log_bias(:)';lab{end+1}=sprintf('%s %s',m{1},c{1}); %#ok<AGROW>
end,end
bar(B');set(gca,'XTickLabel',{'layer 1','layer 2','layer 3','layer 4'});yline(0,'k-');grid on;ylabel('ln(Q_{est}/Q_{true})');legend(lab,'Location','southwest','NumColumns',2);
title('Calibrated layer Q bias (ghost: C1 -0.20@10ms, C2 -0.12@16+0.03@32ms, C3 original)');
exportgraphics(f,fullfile(root,'analysis','confirm_summary.png'),'Resolution',160);close(f);
end
