function make_figure()
root=fileparts(mfilename('fullpath'));r=[jsondecode(fileread(fullfile(root,'analysis','trueq_metrics.json')));jsondecode(fileread(fullfile(root,'analysis','coarse_true_q_metrics.json')))];
dicts={'Cal','CoarseTrueQ','TrueQ','TrueQW'};labels={'calibrated 4-layer Q','true 4-layer mean Q (post hoc)','true per-trace Q','true Q + true wavelet'};
data={'noisy_Moderate','noisy_Strong','clean_primary'};titles={'Moderate noise','Strong noise','noise-free primaries'};
f=figure('Visible','off','Position',[0 0 1500 480]);tiledlayout(1,3,'TileSpacing','compact');
for k=1:3
 nexttile;M=zeros(4,2);
 for i=1:4,for m=1:2,meth={'FullCalRefit','Standard'};s=strcmp({r.dictionary},dicts{i})&strcmp({r.data},data{k})&strcmp({r.method},meth{m});M(i,m)=r(s).overall;end,end
 b=bar(M);b(1).FaceColor=[.6 .6 .6];b(2).FaceColor=[.2 .45 .8];set(gca,'XTickLabel',labels,'XTickLabelRotation',25);ylim([0 .75]);grid on;
 title(titles{k});if k==1,ylabel('overall relative error');legend({'FullCalRefit','Standard (coherence)'},'Location','northeast');end
 if k==3,yline(.30,'r--','registered limit 0.30');end
end
sgtitle('Complex model V5: how much error is left when Q / wavelet are known');
exportgraphics(f,fullfile(root,'analysis','trueq_breakdown.png'),'Resolution',160);close(f);
end
