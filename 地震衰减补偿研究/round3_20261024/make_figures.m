function make_figures()
% Round 3 gathers: shared colour limit per scene (99.5% |target|), no trace normalisation.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'figures');if ~isfolder(out),mkdir(out);end
show={'Target','FullCalRefit','Standard','QScanRefit','Guarded_Standard'};tau=.01;
for scene={'F4_Q-30_Strong','F3_Q+0_Strong','F4_Q+0_Moderate','V5_Q-30_Strong'}
 s=scene{1};p=split(s,'_');r=load(fullfile(root,'results',[s '.mat']),'gathers','scan');g=r.gathers;
 [~,k1]=min(abs(log(r.scan.grid)));gain=(r.scan.objective(k1)-min(r.scan.objective))/r.scan.objective(k1);
 if gain>tau,g.Guarded_Standard=g.QScan_Standard;else,g.Guarded_Standard=g.Standard;end
 T=load(fullfile(root,'dev',[p{1} '_truth.mat']),'tr');tg=T.tr(strcmp({T.tr.name},p{3})).target;
 in=load(fullfile(root,'dev',[p{1} '.mat']),'config');t=in.config.t(:);clim=prctile(abs(tg(:)),99.5);
 f=figure('Visible','off','Position',[0 0 2200 900]);tl=tiledlayout(2,numel(show),'TileSpacing','compact','Padding','compact');
 cmap=[min(1,1+linspace(-1,1,256)'),1-abs(linspace(-1,1,256)'),min(1,1-linspace(-1,1,256)')];
 for row=1:2
  for m=1:numel(show)
   nexttile;
   if row==2&&m==1
    J=r.scan.objective;semilogx(r.scan.grid,(J-min(J))/J(k1),'o-');hold on;xline(1,'k:');xline(r.scan.factor,'r-');
    yline(tau,'b--');xlabel('Q correction factor');ylabel('(J - min J) / J(1)');title('QScan objective','FontSize',9);grid on;continue;
   end
   if strcmp(show{m},'Target'),A=tg;lab='Target';else,A=g.(show{m});lab=show{m};if row==2,A=A-tg;end,end
   imagesc(1:size(tg,2),t,A,[-clim clim]);colormap(gca,cmap);
   if row==2,title(sprintf('%s - target (rel %.3f)',lab,norm(A,'fro')/norm(tg,'fro')),'Interpreter','none','FontSize',9);
   else,title(lab,'Interpreter','none','FontSize',9);end
   xlabel('trace');if m==1,ylabel('time (s)');end
  end
 end
 title(tl,sprintf('%s   (QScan factor %.2f, objective gain %.3g, guard %s)',strrep(s,'_',' '),r.scan.factor,gain,ternary(gain>tau,'applied','not applied')));
 exportgraphics(f,fullfile(out,[s '.png']),'Resolution',170);close(f);
end
end

function r=ternary(c,a,b)
if c,r=a;else,r=b;end
end
