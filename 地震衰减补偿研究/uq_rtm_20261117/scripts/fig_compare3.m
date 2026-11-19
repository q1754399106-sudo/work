function fig_compare3(name,m)
% Descriptive: reference | attenuated (uncompensated, m30) | route-A compensated (UQW1 at Q multiplier m),
% images (top) and image traces every 20 m (bottom). Attenuated and compensated images are scaled to the
% reference energy in 700-1550 m for display. figures/compare3_<name>_m<m>.png
if nargin<1,name='C2';end
if nargin<2,m=1;end
P=rtm_setup();S=load(fullfile(P.cache,['images2_' name '.mat']),'ref','m30','meta');rows=S.meta.combo.m30;
e1=find(strcmp(rows(:,1),'E1'));mm=[rows{e1,2}];z=P.z(P.in{1});x=P.x(P.in{2});iz=z>=700&z<=1550;
ref=prep(S.ref.R0,P);nrm=@(a) a*(norm(ref(iz,:),'fro')/max(norm(a(iz,:),'fro'),eps));
att=nrm(prep(S.m30.R0,P));uq=nrm(prep(S.m30.UQW1(:,e1(mm==m)),P));
pan={'参照（无衰减、无噪声）',ref;'衰减后（未补偿，道上 −30 dB）',att;sprintf('路线 A 补偿后（UQ，m = %g）',m),uq};
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
c=prctile(abs(ref(iz,:)),99,'all');xs=600:20:1800;[~,jx]=min(abs(x(:)-xs),[],1);sc=prctile(abs(ref(iz,jx)),99.5,'all');
f=figure('Visible','off','Position',[0 0 1800 1000]);tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for p=1:3
 nexttile(p);imagesc(x,z(iz),pan{p,2}(iz,:),[-c c]);colormap(gray);axis image;title(pan{p,1});if p==1,ylabel('深度（m）');end
end
for p=1:3
 nexttile(3+p);hold on;A=pan{p,2};zz=z(iz)';
 for j=1:numel(jx)
  s=max(min(A(iz,jx(j))/sc*14,30),-30);x0=xs(j);
  fill([x0+max(s,0);x0*ones(size(s))],[zz;flipud(zz)],'k','EdgeColor','none');plot(x0+s,zz,'k','LineWidth',.3);
 end
 set(gca,'YDir','reverse');ylim([700 1550]);xlim([560 1840]);box on;grid on;xlabel('x（m）');if p==1,ylabel('深度（m）');end
end
exportgraphics(f,fullfile(P.root,'figures',sprintf('compare3_%s_m%03d.png',name,round(100*m))),'Resolution',130);close(f);
end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=10;h=exp(-.5*((-30:30)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),30,1);y;repmat(y(end,:),30,1)],'valid');
y(P.z(P.in{1})<140,:)=0;
end
