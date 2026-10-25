function make_figures()
% Sealed-validation gathers: shared colour limit per scene, no trace normalisation.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'figures');if ~isfolder(out),mkdir(out);end
show={'Target','SoftGroupRefit','v1','Standard','Robust'};
for scene={'V4_Q+0_Strong','V4_Q-30_Strong','V5_Q+0_Strong','V5_Q-30_Strong'}
 s=scene{1};p=split(s,'_');r=load(fullfile(root,'results',[s '.mat']),'gathers');
 T=load(fullfile(root,'truth',[p{1} '_truth.mat']),'tr');tg=T.tr(strcmp({T.tr.name},p{3})).target;
 in=load(fullfile(root,'inputs','sealed',[p{1} '.mat']),'config');t=in.config.t(:);clim=prctile(abs(tg(:)),99.5);
 f=figure('Visible','off','Position',[0 0 2200 900]);tl=tiledlayout(2,numel(show),'TileSpacing','compact','Padding','compact');
 for row=1:2
  for m=1:numel(show)
   nexttile;
   if strcmp(show{m},'Target'),A=tg;lab='Target';if row==2,A=r.gathers.FullCalRefit-tg;lab='FullCalRefit - target';end
   else,A=r.gathers.(show{m});lab=show{m};if row==2,A=A-tg;lab=[lab ' - target'];end,end
   imagesc(1:99,t,A,[-clim clim]);colormap(gca,[min(1,1+linspace(-1,1,256)'),1-abs(linspace(-1,1,256)'),min(1,1-linspace(-1,1,256)')]);
   if row==2,title(sprintf('%s (rel %.3f)',lab,norm(A,'fro')/norm(tg,'fro')),'Interpreter','none','FontSize',9);else,title(lab,'Interpreter','none','FontSize',9);end
   xlabel('trace');if m==1,ylabel('time (s)');end
  end
 end
 title(tl,strrep(s,'_',' '));exportgraphics(f,fullfile(out,[s '.png']),'Resolution',180);close(f);
end
end
