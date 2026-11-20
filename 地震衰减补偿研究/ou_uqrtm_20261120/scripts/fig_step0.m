function fig_step0()
% Step 0 figure: local-field sigma_eff maps (near-offset shot per column) for D2 and D4, and the
% Monte-Carlo check. Writes figures/step0_sigma_eff.png and analysis/step0_sigma_eff.json.
P=ou_setup();maxNumCompThreads(P.threads);ou=P.ou;
names={'D2','D4'};zi=P.z(P.in{1});xi=P.x(P.in{2});xs=P.x(P.ixs);
fig=figure('Visible','off','Position',[50 50 1500 820]);J=struct();
for n=1:2
 M=make_model(names{n},P);Fl=ou_flags(M,P);sloc=ou.s_bg+(ou.s_flag-ou.s_bg)*Fl;
 R=ou_sigma(P,M.Qm,M.vmig,sloc,ou);
 % near-offset value: for every x, the shot closest to x
 [~,kx]=min(abs(xi(:)-xs),[],2);Sg=zeros(numel(zi),numel(xi));
 for j=1:numel(xi),Sg(:,j)=sqrt(R.L2((j-1)*numel(zi)+(1:numel(zi)),kx(j)));end
 % median over shots
 Sm=reshape(sqrt(median(R.L2,2)),numel(zi),numel(xi));
 subplot(2,3,3*(n-1)+1);imagesc(xi,zi,Fl(P.in{1},P.in{2}));axis image;title([names{n} ' 标记区 (s_{loc}=0.6)']);xlabel('x (m)');ylabel('z (m)');
 subplot(2,3,3*(n-1)+2);imagesc(xi,zi,Sg,[0 .3]);axis image;colorbar;title([names{n} ' 局部场 \sigma_{eff}（近炮检距）']);xlabel('x (m)');
 subplot(2,3,3*(n-1)+3);
 wd=zi>=820&zi<=1550;iw=xi>=300&xi<=2100;
 plot(zi,median(Sg(:,iw),2),'k',zi,prctile(Sg(:,iw),95,2),'r--',zi,sqrt(2*ou.Lz./zi)*ou.s_bg,'b:','LineWidth',1.2);hold on;
 yline(.447,'g','标量 UQ 0.447');ylim([0 .5]);xlim([140 1600]);legend('中位','95%','(2L_z/z)^{1/2} s_{bg}','Location','northeast');
 xlabel('z (m)');ylabel('\sigma_{eff}');title([names{n} ' 随深度']);grid on
 J.(names{n})=struct('deep_median',median(Sg(wd,iw),'all'),'deep_p95',prctile(Sg(wd,iw),95,'all'),'deep_max',max(Sg(wd,iw),[],'all'), ...
  'deep_median_allshots',median(Sm(wd,iw),'all'),'deep_max_allshots',max(Sm(wd,iw),[],'all'));
end
exportgraphics(fig,fullfile(P.root,'figures','step0_sigma_eff.png'),'Resolution',110);close(fig);
fid=fopen(fullfile(P.root,'analysis','step0_sigma_eff.json'),'w');fprintf(fid,'%s',jsonencode(J,'PrettyPrint',true));fclose(fid);disp(J.D2);disp(J.D4);
end
