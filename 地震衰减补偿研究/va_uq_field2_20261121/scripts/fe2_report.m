function fe2_report()
% Tables and figures for RESULTS_ZH.md (written after the run; descriptive, no judgement):
% analysis/fe2_tables.md, analysis/fe2_checks.json, figures/fe2_ties.png, figures/fe2_images.png.
[C,P]=fe2_cfg();A=fullfile(C.root,'analysis');R=jsondecode(fileread(fullfile(A,'fe2_eval.json')));
B=jsondecode(fileread(fullfile(A,'fe2_boot.json')));S4=jsondecode(fileread(fullfile(C.s4,'analysis','stage4_field.json')));
S25=jsondecode(fileread(fullfile(C.s4,'analysis','stage4_field_fref25.json')));SA=jsondecode(fileread(fullfile(C.s4,'analysis','stage4_field_amp.json')));
ms=R.m(:)';I=load(fullfile(C.cache,'images_fe2.mat'));I4=load(fullfile(C.s4cache,'images.mat'),'UQW1','R0');I25=load(fullfile(C.s4cache,'images_fref25.mat'),'UQW1');
rel=@(a,b) norm(double(a(:))-double(b(:)))/norm(double(b(:)));
K.R0_vs_stage4_image=rel(I.R0,I4.R0);K.UQW1_100_B1_vs_stage4_image=rel(I.UQW1_100_B1,I4.UQW1);K.UQW1_25_B1_vs_fref25_image=rel(I.UQW1_25_B1,I25.UQW1);
K.R0_tie=R.R0.tie_F5;K.R0_tie_stage4=S4.R0.tie_F5;
K.UQW1_100_B1_tie_maxdiff=max(abs(R.UQW1_100_B1.tie_F5-S4.UQW1.tie_F5));K.UQW1_25_B1_tie_maxdiff=max(abs(R.UQW1_25_B1.tie_F5-S25.UQW1.tie_F5));
% noise estimate: reservoir-window median snr, B300 vs B1
z=P.z(P.in{1});x=P.x(P.in{2});wr=(z>=2600&z<=3200)'&(x>=-2500&x<=1500);
N3=load(fullfile(C.cache,'inv300.mat'),'iv');N1=load(fullfile(C.s4cache,'inv.mat'),'iv');
s3=1./median(N3.iv(wr(:),:),1);s1=1./median(N1.iv(wr(:),:),1);K.f=P.f;K.snr_dB_B300=10*log10(s3);K.snr_dB_B1=10*log10(s1);
fid=fopen(fullfile(A,'fe2_checks.json'),'w');fprintf(fid,'%s',jsonencode(K,'PrettyPrint',true));fclose(fid);
% table
L={};L{end+1}=sprintf('## F-5 井旁相关（R0 带内箱形，fld_eval.m），R0 = %.3f [%.2f, %.2f]\n',R.R0.tie_F5,B.R0.ci90);
L{end+1}=['| 方法 | ' strjoin(arrayfun(@(m) sprintf('m = %g',m),ms,'uni',0),' | ') ' |'];L{end+1}=['|---' repmat('|---',1,numel(ms)) '|'];
rows={'VA','VA（B300）';'UQW1_100','UQW1(100)（B300）';'UQW1_25','UQW1(25)（B300）';'VA_B1','VA（B1，描述）';'UQW1_100_B1','UQW1(100)（B1，复算）';'UQW1_25_B1','UQW1(25)（B1，复算）'};
for r=1:size(rows,1),v=R.(rows{r,1}).tie_F5;L{end+1}=['| ' rows{r,2} ' | ' strjoin(arrayfun(@(a) sprintf('%.3f',a),v,'uni',0),' | ') ' |'];end %#ok<AGROW>
L{end+1}=['| 阶段 4 UQW1(100)（B1，缓存） | ' strjoin(arrayfun(@(a) sprintf('%.3f',a),S4.UQW1.tie_F5,'uni',0),' | ') ' |'];
L{end+1}=['| 阶段 4 UQW1(25)（B1，缓存） | ' strjoin(arrayfun(@(a) sprintf('%.3f',a),S25.UQW1.tie_F5,'uni',0),' | ') ' |'];
L{end+1}=['| 阶段 4 UQW1 只补振幅（B1，缓存） | ' strjoin(arrayfun(@(a) sprintf('%.3f',a),SA.UQW1.tie_F5,'uni',0),' | ') ' |'];
L{end+1}=sprintf('\n## 90%% 自助区间（相关 / 方法 − R0），块长 10 样点，2000 次\n');
L{end+1}=['| 方法 | ' strjoin(arrayfun(@(m) sprintf('m = %g',m),ms,'uni',0),' | ') ' |'];L{end+1}=['|---' repmat('|---',1,numel(ms)) '|'];
for r=1:3
 b=B.(rows{r,1});c=arrayfun(@(k) sprintf('[%.2f, %.2f]',b.ci90(k,1),b.ci90(k,2)),1:numel(ms),'uni',0);d=arrayfun(@(k) sprintf('[%+.2f, %+.2f]',b.diff_ci90(k,1),b.diff_ci90(k,2)),1:numel(ms),'uni',0);
 L{end+1}=['| ' rows{r,2} ' 相关 | ' strjoin(c,' | ') ' |'];L{end+1}=['| ' rows{r,2} ' − R0 | ' strjoin(d,' | ') ' |']; %#ok<AGROW>
end
L{end+1}=sprintf('\n## 储层窗谱质心（周/km；R0 %.2f）\n',R.R0.cen);
L{end+1}=['| 方法 | ' strjoin(arrayfun(@(m) sprintf('m = %g',m),ms,'uni',0),' | ') ' |'];L{end+1}=['|---' repmat('|---',1,numel(ms)) '|'];
for r=1:size(rows,1),v=R.(rows{r,1}).cen;L{end+1}=['| ' rows{r,2} ' | ' strjoin(arrayfun(@(a) sprintf('%.2f',a),v,'uni',0),' | ') ' |'];end %#ok<AGROW>
fid=fopen(fullfile(A,'fe2_tables.md'),'w');fprintf(fid,'%s\n',L{:});fclose(fid);
% figures
set(0,'DefaultAxesFontName','Microsoft YaHei','DefaultTextFontName','Microsoft YaHei','DefaultAxesFontSize',9);
f=figure('Visible','off','Position',[0 0 1400 450]);tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on;cl=lines(6);
yregion(0.60,0.66,'FaceColor',[.85 .95 .85],'EdgeColor','none','DisplayName','P4 区间 0.60–0.66');
yline(R.R0.tie_F5,'k--','R0','HandleVisibility','off');
for r=1:3,b=B.(rows{r,1});errorbar(ms,R.(rows{r,1}).tie_F5,R.(rows{r,1}).tie_F5(:)-b.ci90(:,1),b.ci90(:,2)-R.(rows{r,1}).tie_F5(:),'-o','Color',cl(r,:),'LineWidth',1.2+1.2*(r==1),'DisplayName',rows{r,2});end
plot(ms,S4.UQW1.tie_F5,':s','Color',cl(2,:),'DisplayName','阶段 4 UQW1(100) B1');plot(ms,S25.UQW1.tie_F5,':s','Color',cl(3,:),'DisplayName','阶段 4 UQW1(25) B1');
set(gca,'XScale','log');grid on;xlabel('Q 乘子 m（名义 Q = 120）');ylabel('F-5 井旁相关');title('F-5 箱形带相关（误差棒：90% 块自助）');legend('Location','southeast');
nexttile;plot(P.f,K.snr_dB_B1,'k',P.f,K.snr_dB_B300,'r','LineWidth',1.3);grid on;xlabel('频率（Hz）');ylabel('dB');legend('B = 1（阶段 4）','B = 300（本次）');title('储层窗成像域信噪比估计（中位数）');
exportgraphics(f,fullfile(C.root,'figures','fe2_ties.png'),'Resolution',130);close(f);
pr=@(v) prep(v,P);zz=z>=1500&z<=3400;ref=pr(I.R0);c=prctile(abs(ref(zz,:)),99,'all');
pans={'R0',1,'不补偿 R0';'VA',2,'VA，m = 0.6';'VA',5,'VA，m = 1';'UQW1_25',2,'UQW1(25)，m = 0.6';'UQW1_25',5,'UQW1(25)，m = 1';'UQW1_100',5,'UQW1(100)，m = 1'};
f=figure('Visible','off','Position',[0 0 1800 700]);tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for p=1:size(pans,1)
 a=I.(pans{p,1});a=pr(a(:,min(pans{p,2},size(a,2))));a=a*(norm(ref(zz,:),'fro')/max(norm(a(zz,:),'fro'),eps));
 nexttile;imagesc(x/1000,z(zz),a(zz,:),[-c c]);colormap(gray);hold on;plot([-1.275 -1.118],[2570 3170],'r-','LineWidth',1.5);
 Aa=load(fullfile(C.cache,'anchor.mat'),'za');plot(x/1000,Aa.za,'c-');title(pans{p,3});
end
exportgraphics(f,fullfile(C.root,'figures','fe2_images.png'),'Resolution',110);close(f);
disp(K)
end

function y=prep(img,P)
y=reshape(double(img),P.nzi,P.nxi);s=4;h=exp(-.5*((-12:12)/s).^2);h=h/sum(h);
y=y-conv2(h',1,[repmat(y(1,:),12,1);y;repmat(y(end,:),12,1)],'valid');y(P.z(P.in{1})<150,:)=0;
end
