function make_figures()
% Shared colour limits per scene (99.5% |amplitude| of target), no trace normalisation.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'figures');if ~isfolder(out),mkdir(out);end
show={'Target','SoftGroupRefit','CoherenceSoftRefit','PrunedFullCalRefit'};
for scene={'V1_Q+0_Strong','V2_Q+0_Strong','V2_Q-30_Strong','V3_Q+0_Strong','V3_Q+0_Moderate'}
 s=scene{1};parts=split(s,'_');r=load(fullfile(root,'results',[s '.mat']),'gathers','diagnostics');
 T=load(fullfile(root,'truth',[parts{1} '_truth.mat']),'tr');tg=T.tr(strcmp({T.tr.name},parts{3})).target;
 in=load(fullfile(root,'inputs','validation',[parts{1} '.mat']),'config');t=in.config.t(:);
 clim=prctile(abs(tg(:)),99.5);f=figure('Visible','off','Position',[0 0 1800 900]);
 tl=tiledlayout(2,numel(show),'TileSpacing','compact','Padding','compact');
 for row=1:2
  for m=1:numel(show)
   nexttile;
   if strcmp(show{m},'Target'),A=tg;else,A=r.gathers.(show{m});end
   if row==2
    if strcmp(show{m},'Target'),A=r.gathers.FullCalRefit-tg;ttl='FullCalRefit - target';else,A=A-tg;ttl=[show{m} ' - target'];end
   else,ttl=show{m};end
   imagesc(1:99,t,A,[-clim clim]);colormap(gca,gray_red());
   sat=mean(abs(A(:))>clim)*100;
   if row==2,e=norm(A,'fro')/norm(tg,'fro');title(sprintf('%s  (rel %.3f, clip %.1f%%)',ttl,e,sat),'Interpreter','none','FontSize',9);
   else,title(sprintf('%s  (clip %.1f%%)',ttl,sat),'Interpreter','none','FontSize',9);end
   xlabel('trace');if m==1,ylabel('time (s)');end
  end
 end
 title(tl,strrep(s,'_',' '));exportgraphics(f,fullfile(out,[s '.png']),'Resolution',200);close(f);
end
% Coherence map for one scene: c of each group placed at its centre sample.
r=load(fullfile(root,'results','V2_Q+0_Strong.mat'),'diagnostics');
f=figure('Visible','off','Position',[0 0 900 700]);hold on;
cc=r.diagnostics.coherence;
for block=1:9
 c=load(fullfile(root,'cache',sprintf('V2_Q+0_Strong_%02d.mat',block)),'initial');
 for q=1:11
  j=(block-1)*11+q;s=c.initial.support{q}(:);if isempty(s),continue;end
  g=cumsum([1;diff(s)>4]);cen=accumarray(g,s,[],@mean);
  scatter(j*ones(size(cen)),(cen-1)*.004,14,cc{j},'filled');
 end
end
set(gca,'YDir','reverse');colormap(parula);colorbar;caxis([0 1]);xlabel('trace');ylabel('time (s)');
title('V2 Strong Q0: neighbour support c per reflection group');box on;
exportgraphics(f,fullfile(out,'V2_coherence_map.png'),'Resolution',200);close(f);
end

function m=gray_red()
x=linspace(-1,1,256)';m=[min(1,1+x),1-abs(x),min(1,1-x)];
end
