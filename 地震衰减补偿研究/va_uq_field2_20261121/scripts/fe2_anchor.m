function A=fe2_anchor()
% Anchor horizon (PLAN_ZH.md section 3), fixed before any score. Rule:
%  1. in every column of the stage-4 VTI 2008 velocity section (P.vp, unsmoothed, 10 m), the 40 m centred
%     vertical velocity jump d(z) = vp(z+20) - vp(z-20);
%  2. z_raw(x) = depth of the largest d in 2150-2400 m (the layer boundary at about 2.3 km above the
%     reservoir; the next boundary, top of the 4.2-4.5 km/s layer at 2.49-2.56 km, is too close to the
%     F-5 tie window 2570-3170 m and the reservoir window 2600-3200 m);
%  3. z_a(x) = 21-point (210 m) running median of z_raw along the line, rounded to the 10 m grid;
%     PML columns take the edge values. 1/Q_c = 0 above z_a(x) (z < z_a), Q_c = m*120 below.
% Also reports how well the uncompensated stage-4 image R0 shows a reflection on z_a (envelope within
% +-20 m against the median envelope of 1900-2500 m in the same column), the velocity step, and the
% clearance to the F-5 tie window and the reservoir window. Writes cache/anchor.mat,
% analysis/anchor.json, figures/anchor.png.
[C,P]=fe2_cfg();z=P.z(P.in{1});x=P.x(P.in{2});nx=numel(x);vp=P.vp(P.in{1},P.in{2});
h=round(C.anchor.jump_half/P.dx);d=nan(size(vp));d(1+h:end-h,:)=vp(1+2*h:end,:)-vp(1:end-2*h,:);
w=find(z>=C.anchor.zwin(1)&z<=C.anchor.zwin(2));zr=zeros(1,nx);jmp=zr;
for c=1:nx,[jmp(c),k]=max(d(w,c));zr(c)=z(w(k));end
za=round(movmedian(zr,C.anchor.medlen,'Endpoints','shrink')/P.dx)*P.dx;
za_full=[za(1)*ones(1,P.npml) za za(end)*ones(1,P.nx-P.npml-nx)];
iz=round((za-z(1))/P.dx)+1;va=zeros(1,nx);vb=va;
for c=1:nx,va(c)=vp(max(iz(c)-3,1),c);vb(c)=vp(min(iz(c)+3,numel(z)),c);end
% reflection strength of R0 on the horizon
I=load(fullfile(C.s4cache,'images.mat'),'R0');R0=prep(I.R0,P);e=abs(hilbert(R0));er=zeros(1,nx);
kb=z>=1900&z<=2500;for c=1:nx,er(c)=max(e(abs(z-za(c))<=20,c))/median(e(kb,c));end
wx=x>=-2500&x<=1500;wf=x>=-1280&x<=-1110;                  % evaluation window, F-5 tie columns
A=struct('rule','max 40 m vp jump in 2150-2400 m per column, 21-point running median, 10 m grid', ...
 'za_line',[min(za) max(za)],'za_evalwin',[min(za(wx)) max(za(wx))],'za_F5',[min(za(wf)) max(za(wf))], ...
 'jump_median',median(jmp(wx)),'jump_min',min(jmp(wx)),'v_above_median',median(va(wx)),'v_below_median',median(vb(wx)), ...
 'env_ratio_median',median(er(wx)),'env_ratio_p10',prctile(er(wx),10),'clear_to_F5_tie_top',2570-max(za(wf)), ...
 'clear_to_reservoir_top',2600-max(za(wx)),'raw_minus_median_absmax',max(abs(zr-za)));
save(fullfile(C.cache,'anchor.mat'),'za','za_full','zr','jmp','x');
fid=fopen(fullfile(C.root,'analysis','anchor.json'),'w');fprintf(fid,'%s',jsonencode(A,'PrettyPrint',true));fclose(fid);
disp(A)
% figure
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
f=figure('Visible','off','Position',[0 0 1500 520]);tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile([1 2]);zz=z>=1500&z<=3400;c=prctile(abs(R0(zz,:)),99,'all');imagesc(x/1000,z(zz),R0(zz,:),[-c c]);colormap(gray);hold on
plot(x/1000,za,'c-','LineWidth',1.2);plot([-1.275 -1.118],[2570 3170],'r-','LineWidth',1.5);
rectangle('Position',[-2.5 2600 4 600],'EdgeColor','y','LineStyle','--');
xlabel('测线坐标（km，平台为 0）');ylabel('深度（m）');title('阶段 4 不补偿像 R0、锚点 z_a(x)（青）、F-5（红）、储层窗（黄）');
nexttile;[~,jw]=min(abs(x+1200));plot(vp(:,jw),z,'k');hold on;plot(P.vmig(P.in{1},P.in{2}(jw)),z,'b');yline(za(jw),'c-','z_a');
set(gca,'YDir','reverse');ylim([1500 3400]);grid on;xlabel('速度（m/s）');legend('VTI 2008','偏移速度（平滑）','Location','southwest');title(sprintf('x = %.2f km（F-5 附近）',x(jw)/1000));
exportgraphics(f,fullfile(C.root,'figures','anchor.png'),'Resolution',120);close(f);
end

function y=prep(img,P)
% as fld_eval.m (40 m smoothing removed, z < 150 m zeroed)
y=reshape(double(img),P.nzi,P.nxi);s=4;h=exp(-.5*((-12:12)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),12,1);y;repmat(y(end,:),12,1)],'valid');y(P.z(P.in{1})<150,:)=0;
end
