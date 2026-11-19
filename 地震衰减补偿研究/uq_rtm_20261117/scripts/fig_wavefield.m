function fig_wavefield(name)
% Descriptive: pressure snapshots of one shot (x = 1200 m) of model name, elastic (Q = Inf) against
% viscoacoustic (GSLS), at 0.3 / 0.6 / 0.9 s, plus the model. figures/wavefield_<name>.png
if nargin<1,name='D2';end
P=rtm_setup();maxNumCompThreads(4);M=make_model(name,P);
ts=[.3 .6 .9];nt=round(max(ts)/P.fd.dt)+10;steps=round(ts/P.fd.dt);
t=(1:nt)'*P.fd.dt-P.fd.t0;a=(pi*P.fd.fpeak*t).^2;q=(1-2*a).*exp(-a);
opt=struct('dx',P.dx,'dt',P.fd.dt,'nt',nt,'fref',P.fd.fref,'npml',P.npml,'src',[P.izs P.npml+1+300],'rec',[P.izs P.npml+1+300],'snap',steps);
[~,Sa]=fd2d_visco_snap(M,q,opt);Me=M;Me.Q(:)=Inf;[~,Se]=fd2d_visco_snap(Me,q,opt);
iz=P.in{1};ix=P.in{2};z=P.z(iz);x=P.x(ix);
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
f=figure('Visible','off','Position',[0 0 1500 900]);tl=tiledlayout(3,4,'TileSpacing','compact','Padding','loose');
nexttile;imagesc(x,z,M.v(iz,ix));axis image;colorbar;title('速度（m/s）');ylabel('深度（m）');
nexttile(5);imagesc(x,z,M.Q(iz,ix),[15 120]);axis image;colorbar;title('Q');ylabel('深度（m）');
nexttile(9);axis off;text(0,.6,{'压力波场快照：','无衰减 / 有衰减 / 两者之差','',sprintf('模型 %s，炮点 x = 1200 m',name),'Ricker 主频 25 Hz','二维黏声有限差分（GSLS）','每张快照用无衰减快照的','99% 分位数统一定标'},'FontSize',10);
for k=1:3
 e=Se(iz,ix,k);g=Sa(iz,ix,k);c=prctile(abs(e(:)),99);
 nexttile(4*(k-1)+2);imagesc(x,z,e,[-c c]);axis image;title(sprintf('无衰减，t = %.1f s',ts(k)));
 nexttile(4*(k-1)+3);imagesc(x,z,g,[-c c]);axis image;title(sprintf('有衰减，t = %.1f s',ts(k)));
 nexttile(4*(k-1)+4);imagesc(x,z,e-g,[-c c]);axis image;title(sprintf('差（衰减掉的部分），t = %.1f s',ts(k)));
end
colormap(gca,'gray');for ax=findall(f,'Type','axes')',if ax.Title.String~="速度（m/s）"&&ax.Title.String~="Q",colormap(ax,gray);end,end
exportgraphics(f,fullfile(P.root,'figures',['wavefield_' name '.png']),'Resolution',140);close(f);
end
