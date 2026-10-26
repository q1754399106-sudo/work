function plot_research_models(resume)
% MATLAB R2024b. All scientific rendering is performed by MATLAB.
% Run: addpath('.../model_visualization_20261026'); plot_research_models
root=fileparts(mfilename('fullpath')); inp=fullfile(root,'inputs'); out=fullfile(root,'figures');
if nargin<1,resume=false;end
setappdata(groot,'model_visualization_resume',resume);
if ~isfolder(out),mkdir(out);end
diary(fullfile(root,'MATLAB_run.log')); cleanup=onCleanup(@()diary('off'));
maxNumCompThreads(4);
set(groot,'defaultFigureVisible','off','defaultAxesFontName','Microsoft YaHei', ...
 'defaultTextFontName','Microsoft YaHei','defaultAxesFontSize',10, ...
 'defaultTextInterpreter','none','defaultAxesTickLabelInterpreter','none', ...
 'defaultLegendInterpreter','none');
runtime=fullfile(root,'cache','runtime');if ~isfolder(runtime),mkdir(runtime);end
% Extract geometry-only code verbatim. Do not rerun noise, calibration or wave solvers.
s=fileread(fullfile(inp,'calq_confirm_20261025','core','gen_complex_model.m'));
s=extractBefore(s,'N=4096;');s=strrep(s,'function model=gen_complex_model(g)','function model=geometry_only(g)');
s=[s newline 'model=struct(''x'',x,''z'',z,''vp'',vp,''rho'',rho,''Q'',Q,''strata'',strata,''wedge'',wedge,''channel'',channel,''gas'',gas);' newline 'end'];
writeutf(fullfile(runtime,'geometry_only.m'),s);
s=fileread(fullfile(inp,'fd2d_20261025','scripts','fd_setup.m'));
a=strfind(s,'P.root=');b=strfind(s,'% Sweep constants:');s=[s(1:a-1) s(b:end)];
s=strrep(s,'function P=fd_setup()','function P=fd_display_setup()');
writeutf(fullfile(runtime,'fd_display_setup.m'),s);
copyfile(fullfile(inp,'fd2d_20261025','scripts','fd_geology.m'),fullfile(runtime,'fd_geology.m'));
addpath(runtime,fullfile(inp,'calq_confirm_20261025','source'));
records=struct('name',{},'family',{},'png',{},'caption',{});
sourceData=struct; checks=struct('matlab',version,'complex_V3_max_error',NaN,'H_base_match',false(1,6));

%% Ten four-layer configurations, including the original control.
spec=jsondecode(fileread(fullfile(inp,'calq_confirm_20261025','inputs','p1_experiment_spec.json')));
c=config_model();c.name='OriginalControl';c.split='control';models={c};
splits={'historical','development','heldout_new_geology'};
for k=1:numel(splits)
 ss=spec.models.(splits{k});for j=1:numel(ss),cc=ss(j);cc.split=splits{k};cc.h=diff([0 reshape(cc.event_t,1,[])]).*reshape(cc.vp(1:3),1,[])/2;models{end+1}=cc;end %#ok<AGROW>
end
cols=lines(4);zmax=1800;
for prop={'vp','Q'}
 p=prop{1};f=newfig([1500 720]);tl=tiledlayout(f,2,5,'TileSpacing','compact','Padding','compact');
 for k=1:numel(models)
  cc=models{k};ax=nexttile(tl);bounds=[0 cumsum(cc.h(:)') zmax];vals=cc.(p)(:);
  for L=1:4,patch(ax,[0 vals(L) vals(L) 0],[bounds(L) bounds(L) bounds(L+1) bounds(L+1)],cols(L,:),'FaceAlpha',.10,'EdgeColor','none');hold(ax,'on');end
  [xx,zz]=stepcoords(vals,bounds);plot(ax,xx,zz,'-','Color',[.08 .20 .32],'LineWidth',1.7);
  for L=2:4,yline(ax,bounds(L),':','Color',[.5 .5 .5]);end
  set(ax,'YDir','reverse');ylim(ax,[0 zmax]);grid(ax,'on');title(ax,strrep(cc.name,'_',' '),'FontSize',10);
  if strcmp(p,'vp'),xlim(ax,[0 4500]);xlabel(ax,'P 波速度 (m/s)');else,xlim(ax,[0 180]);xlabel(ax,'品质因子 Q');end
  ylabel(ax,'深度 (m)');
 end
 title(tl,['四层地质模型 · ' p_label(p)],'FontSize',18,'FontWeight','bold');
 subtitle(tl,'前三层厚度由给定双程时间和层速度换算；第四层为半空间，仅展示至 1800 m','FontSize',10);
 records(end+1)=deliver(f,out,['01_layered_' p],'四层模型',sprintf('10 个配置；%s；共同坐标范围；水平线为层界面。',p_label(p))); %#ok<AGROW>
end
sourceData.layered=models;

%% Sparse reflection-time models V1/V2/V4: do not invent depth sections.
for nn={'V1','V2','V4'}
 n=nn{1};[S,T]=loadV(inp,n);R=T.reflectivity;t=S.config.t(:);x=1:size(R,2);
 target=T.tr(1).target;D=S.obs(1).D;lim=max(abs([target(:);D(:)]));
 f=newfig([1250 980]);tl=tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
 ax=nexttile(tl);imageplot(ax,x,t,R,'反射系数（无量纲）',[-max(abs(R(:))) max(abs(R(:)))],divmap());title(ax,'(a) 原始反射系数模型');
 ax=nexttile(tl);imageplot(ax,x,t,target,'相对振幅',[-lim lim],divmap());title(ax,'(b) 无衰减目标道集');
 ax=nexttile(tl);imageplot(ax,x,t,D,'相对振幅',[-lim lim],divmap());title(ax,'(c) 含衰减观测 · Moderate (14 dB)');
 axlist=findall(f,'Type','axes');for a=axlist',xlabel(a,'道号');ylabel(a,'双程时间 (s)');end
 title(tl,[n ' · 稀疏反射模型'],'FontSize',18);subtitle(tl,'反射系数在时间域定义；(b)(c) 共用对称色标，无逐道归一化、无自动增益','FontSize',10);
 records(end+1)=deliver(f,out,['02_' n '_reflectivity'],'稀疏反射模型','时间域反射系数及配套道集；保留原始极性和幅度。'); %#ok<AGROW>
 sourceData.(n)=struct('reflectivity',R,'t',t,'Q',S.config.Q,'vp',S.config.vp,'h',S.config.h);
end

%% Complex physical sections V3/V5 from archived geometry and original RNG.
complexModels=cell(1,2);
for k=1:2
 ns={'V3','V5'};n=ns{k};[S,T]=loadV(inp,n);G=geometry_only(T.g);
 if isfield(T,'mt'),checks.complex_V3_max_error=max([max(abs(G.vp-T.mt.vp),[],'all'),max(abs(G.Q-T.mt.Q),[],'all')]);assert(checks.complex_V3_max_error==0);end
 complexModels{k}=G;sourceData.(n)=G;
 f=newfig([1350 930]);tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
 ax=nexttile(tl);imageplot(ax,G.x,G.z,G.vp,'P 波速度 (m/s)',[1750 3650],parula(256));feature_contours(ax,G);title(ax,'(a) 速度与地质构造');xlabel(ax,'横向距离 (m)');ylabel(ax,'深度 (m)');
 ax=nexttile(tl);imageplot(ax,G.x,G.z,G.Q,'品质因子 Q',[18 160],parula(256));feature_contours(ax,G);title(ax,'(b) Q 与低 Q 异常体');xlabel(ax,'横向距离 (m)');ylabel(ax,'深度 (m)');
 target=T.tr(1).target;D=S.obs(1).D;lim=max(abs([target(:);D(:)]));
 ax=nexttile(tl);imageplot(ax,G.x,S.config.t,target,'相对振幅',[-lim lim],divmap());title(ax,'(c) 无衰减目标道集');xlabel(ax,'横向距离 (m)');ylabel(ax,'双程时间 (s)');
 ax=nexttile(tl);imageplot(ax,G.x,S.config.t,D,'相对振幅',[-lim lim],divmap());title(ax,'(d) 含衰减、噪声及干扰的观测');xlabel(ax,'横向距离 (m)');ylabel(ax,'双程时间 (s)');
 title(tl,[n ' · 复杂细层地质模型'],'FontSize',18);
 subtitle(tl,'白实线：河道；黑虚线：楔状体；白点线：含气异常。逐列正演，附加运动学绕射/多次波；非二维全波场','FontSize',10);
 records(end+1)=deliver(f,out,['03_' n '_complex'],'复杂地质模型','速度和 Q 剖面；河道、楔状体、含气异常轮廓；目标和观测共用幅度色标。'); %#ok<AGROW>
end

%% Kennett F1-F4 and U1-U4: exact archived logs, original lateral rule.
fineModels=cell(1,8);fnames=[arrayfun(@(k)sprintf('F%d',k),1:4,'UniformOutput',false),arrayfun(@(k)sprintf('U%d',k),1:4,'UniformOutput',false)];
for k=1:8
 n=fnames{k};
 if n(1)=='F',basepath=fullfile(inp,'round3_20261024','dev');S=load(fullfile(basepath,[n '.mat']));T=load(fullfile(basepath,[n '_truth.mat']));
 else,S=load(fullfile(inp,'b0_world_truth_20261024','inputs',[n '.mat']));T=load(fullfile(inp,'b0_world_truth_20261024','truth',[n '_truth.mat']));end
 z=(0:2:2600)';G=fine_section(T.base,z);fineModels{k}=G;sourceData.(n)=G;
 f=newfig([1350 870]);tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
 ax=nexttile(tl);imageplot(ax,1:99,z,G.vp,'P 波速度 (m/s)',[1700 4400],parula(256));title(ax,'(a) 细层速度');xlabel(ax,'道号');ylabel(ax,'深度 (m)');
 ax=nexttile(tl);imageplot(ax,1:99,z,G.Q,'品质因子 Q',[0 260],parula(256));title(ax,'(b) 细层 Q');xlabel(ax,'道号');ylabel(ax,'深度 (m)');
 ax=nexttile(tl);plot(ax,G.vp(:,50),z,'Color',[.05 .35 .65],'LineWidth',1);set(ax,'YDir','reverse');ylim(ax,[0 2600]);grid(ax,'on');xlabel(ax,'P 波速度 (m/s)');ylabel(ax,'深度 (m)');title(ax,'(c) 第 50 道速度测井');
 ax=nexttile(tl);plot(ax,G.Q(:,50),z,'Color',[.65 .25 .05],'LineWidth',1);set(ax,'YDir','reverse');ylim(ax,[0 2600]);grid(ax,'on');xlabel(ax,'品质因子 Q');ylabel(ax,'深度 (m)');title(ax,'(d) 第 50 道 Q 测井');
 title(tl,[n ' · Kennett 细层模型'],'FontSize',18);subtitle(tl,'原始细层日志 + 原横向厚度/尖灭规则；按 2 m 深度网格取所属层值；无平滑；横坐标无米制定义','FontSize',10);
 records(end+1)=deliver(f,out,['04_' n '_Kennett'],'Kennett 细层模型','由归档 base 日志按原 lateral 规则重建；2 m 显示采样；第 50 道测井。'); %#ok<AGROW>
end
for prop={'vp','Q'}
 p=prop{1};f=newfig([1450 1500]);tl=tiledlayout(f,4,2,'TileSpacing','compact','Padding','compact');
 for k=1:8,ax=nexttile(tl);if strcmp(p,'vp'),lim=[1700 4400];else,lim=[0 260];end
  imageplot(ax,1:99,fineModels{k}.z,fineModels{k}.(p),p_label(p),lim,parula(256));title(ax,fnames{k});xlabel(ax,'道号');ylabel(ax,'深度 (m)');
 end
 title(tl,['F1–F4 / U1–U4 · ' p_label(p)],'FontSize',18);subtitle(tl,'各面板共用色标范围；原始细层结构，未平滑','FontSize',10);
 records(end+1)=deliver(f,out,['05_Kennett_atlas_' p],'模型总览','八个 Kennett 模型，同变量统一色标。'); %#ok<AGROW>
end

%% True 2-D finite-difference H1-H6 physical model, original Backus grid.
P=fd_display_setup();hmodels=cell(1,6);
maxHQ=0;
for k=1:6
 n=sprintf('H%d',k);G0=fd_geology(P.models(k),P);T=load(fullfile(inp,'fd2d_20261025','truth',[n '_truth.mat']));
 checks.H_base_match(k)=isequaln(G0.base,T.base);assert(checks.H_base_match(k),'H base mismatch');
 ix=G0.x>=0&G0.x<=max(G0.trace_x);iz=G0.z>=0&G0.z<=P.fd.zmax;
 G=struct('x',G0.x(ix),'z',G0.z(iz),'vp',G0.v(iz,ix),'Q',G0.Q(iz,ix),'rho',G0.rho(iz,ix));hmodels{k}=G;sourceData.(n)=G;
 maxHQ=max(maxHQ,max(G.Q,[],'all'));
end
hQlim=[0 ceil(maxHQ/10)*10];checks.H_Q_color_limits=hQlim;
for k=1:6
 n=sprintf('H%d',k);G=hmodels{k};
 f=newfig([1350 730]);tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
 ax=nexttile(tl);imageplot(ax,G.x,G.z,G.vp,'P 波速度 (m/s)',[1700 4400],parula(256));title(ax,'(a) 有效速度');xlabel(ax,'横向距离 (m)');ylabel(ax,'深度 (m)');
 ax=nexttile(tl);imageplot(ax,G.x,G.z,G.Q,'品质因子 Q',hQlim,parula(256));title(ax,'(b) 有效 Q');xlabel(ax,'横向距离 (m)');ylabel(ax,'深度 (m)');
 title(tl,[n ' · 二维黏声有限差分模型'],'FontSize',18);
 subtitle(tl,sprintf('原始 4 m 网格，层内 0.25 m 子采样的 Backus 等效参数；断层落差 %d m；仅展示物理区，PML 未画出',P.models(k).throw),'FontSize',10);
 records(end+1)=deliver(f,out,['06_' n '_FD2D'],'二维有限差分模型','原始 fd_geology 与固定随机种子重建，base 与归档真值逐字段一致；剔除吸收边界。'); %#ok<AGROW>
 fprintf('Rendered %s\n',n);
end
for prop={'vp','Q'}
 p=prop{1};f=newfig([1450 1250]);tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
 for k=1:6,ax=nexttile(tl);if strcmp(p,'vp'),lim=[1700 4400];else,lim=hQlim;end
  imageplot(ax,hmodels{k}.x,hmodels{k}.z,hmodels{k}.(p),p_label(p),lim,parula(256));title(ax,sprintf('H%d',k));xlabel(ax,'横向距离 (m)');ylabel(ax,'深度 (m)');
 end
 title(tl,['H1–H6 · 二维有限差分模型 · ' p_label(p)],'FontSize',18);subtitle(tl,'同变量统一色标；仅展示物理区；坐标为米','FontSize',10);
 records(end+1)=deliver(f,out,['07_FD2D_atlas_' p],'模型总览','六个独立二维模型，同变量统一色标。'); %#ok<AGROW>
end

% One compact overview: four model families, clearly separated physical axes.
f=newfig([1450 950]);tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);cc=models{1};bounds=[0 cumsum(cc.h(:)') 1800];[xx,zz]=stepcoords(cc.vp,bounds);plot(ax,xx,zz,'-','Color',[.05 .3 .55],'LineWidth',2);set(ax,'YDir','reverse');xlim(ax,[1500 3500]);ylim(ax,[0 1800]);grid(ax,'on');xlabel(ax,'P 波速度 (m/s)');ylabel(ax,'深度 (m)');title(ax,'(a) 原始四层模型 · OriginalControl');
ax=nexttile(tl);G=complexModels{2};imageplot(ax,G.x,G.z,G.vp,'P 波速度 (m/s)',[1700 4400],parula(256));feature_contours(ax,G);xlabel(ax,'横向距离 (m)');ylabel(ax,'深度 (m)');title(ax,'(b) 复杂地层 · V5');
ax=nexttile(tl);G=fineModels{3};imageplot(ax,1:99,G.z,G.vp,'P 波速度 (m/s)',[1700 4400],parula(256));xlabel(ax,'道号');ylabel(ax,'深度 (m)');title(ax,'(c) Kennett 细层 · F3');
ax=nexttile(tl);G=hmodels{1};imageplot(ax,G.x,G.z,G.vp,'P 波速度 (m/s)',[1700 4400],parula(256));xlabel(ax,'横向距离 (m)');ylabel(ax,'深度 (m)');title(ax,'(d) 二维有限差分 · H1');
title(tl,'地震衰减补偿研究 · 主要物理模型','FontSize',20,'FontWeight','bold');subtitle(tl,'读取归档模型和参数，MATLAB R2024b 绘制；二维剖面 (b)(c)(d) 共用速度色标','FontSize',11);
records(end+1)=deliver(f,out,'00_research_models_overview','研究模型总览','四个物理模型家族的代表。V1/V2/V4 为独立时间域反射模型，另有专图。');

save(fullfile(root,'model_display_data.mat'),'sourceData','checks','-v7');
writeutf(fullfile(root,'VALIDATION.json'),jsonencode(checks,PrettyPrint=true));
writeutf(fullfile(root,'FIGURE_INDEX.json'),jsonencode(records,PrettyPrint=true));
write_gallery(root,records);
fprintf('COMPLETE: %d figures, 29 model configurations.\n',numel(records));
end

function f=newfig(sz)
f=figure('Color','w','Units','pixels','Position',[40 40 sz],'Visible','off');
end
function label=p_label(p)
if strcmp(p,'vp'),label='P 波速度 (m/s)';else,label='品质因子 Q（无量纲）';end
end
function [x,z]=stepcoords(v,b)
v=v(:)';x=repelem(v,2);z=reshape([b(1:end-1);b(2:end)],1,[]);
end
function imageplot(ax,x,y,A,lab,lim,map)
assert(all(isfinite(A),'all'));assert(min(A,[],'all')>=lim(1)-1e-10&&max(A,[],'all')<=lim(2)+1e-10,'Color range would clip data');
imagesc(ax,x,y,A,lim);set(ax,'YDir','reverse','Layer','top');colormap(ax,map);cb=colorbar(ax);cb.Label.String=lab;cb.Label.Interpreter='none';
end
function map=divmap()
map=interp1([0 .5 1],[.12 .32 .56;1 1 1;.65 .17 .18],linspace(0,1,257));
end
function feature_contours(ax,G)
hold(ax,'on');contour(ax,G.x,G.z,double(G.channel),[.5 .5],'w-','LineWidth',1.1);contour(ax,G.x,G.z,double(G.wedge),[.5 .5],'k--','LineWidth',1.1);contour(ax,G.x,G.z,double(G.gas),[.5 .5],'w:','LineWidth',1.4);
end
function [S,T]=loadV(inp,n)
if any(strcmp(n,{'V1','V2','V3'})),S=load(fullfile(inp,'coherence_soft_20261023','inputs','validation',[n '.mat']));T=load(fullfile(inp,'coherence_soft_20261023','truth',[n '_truth.mat']));
else,S=load(fullfile(inp,'round3_20261024','dev',[n '.mat']));T=load(fullfile(inp,'round3_20261024','dev',[n '_truth.mat']));end
end
function G=fine_section(b,z)
V=zeros(numel(z),99);Q=V;
for j=1:99
 u=(j-1)/98;stretch=[1+.06*u 1+.05*sin(2*pi*u) 1-.05*u 1];hs=b.h(:)'.*stretch(b.macro(:)');
 hs(b.gas)=b.h(b.gas)*min(max((90-j)/50,0),1);keep=hs>0;
 v=[b.cap_v b.v(keep)];q=[b.cap_Q b.Q(keep)];h=[b.cap_h hs(keep)];v(end+1)=v(end)*1.05;q(end+1)=q(end);
 ii=discretize(z,[0 cumsum(h) Inf]);V(:,j)=v(ii);Q(:,j)=q(ii);
end
G=struct('z',z,'trace',1:99,'vp',V,'Q',Q);
end
function rec=deliver(f,out,stem,family,caption)
rec=struct('name',stem,'family',family,'png',['figures/' stem '.png'],'caption',caption);
if getappdata(groot,'model_visualization_resume')&&all(cellfun(@(ext)isfile(fullfile(out,[stem ext])),{'.png','.pdf','.fig'}))
 close(f);fprintf('Verified existing %s\n',stem);return;
end
drawnow;exportgraphics(f,fullfile(out,[stem '.png']),'Resolution',200,'BackgroundColor','white');
exportgraphics(f,fullfile(out,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
savefig(f,fullfile(out,[stem '.fig']),'compact');close(f);
fprintf('Saved %s\n',stem);
end
function writeutf(p,s)
fid=fopen(p,'w','n','UTF-8');assert(fid>=0);fwrite(fid,unicode2native(s,'UTF-8'),'uint8');fclose(fid);
end
function write_gallery(root,rec)
html=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>研究模型图册</title>' ...
 '<style>body{font-family:Microsoft YaHei,Arial;background:#edf1f5;color:#203046;margin:28px auto;max-width:1350px;padding:0 18px}h1{font-size:28px}section{background:white;padding:18px;margin:22px 0;border-radius:12px}img{width:100%;height:auto}a{color:#165a99;margin-right:18px}p{line-height:1.7}</style>' ...
 '<h1>地震衰减补偿研究 · MATLAB 模型图册</h1><p>29 个主要模型配置：四层模型 10 个，V1–V5，F1–F4，U1–U4，H1–H6。每张图均提供 PNG、PDF 和可编辑 FIG。原始数据及校验见 SOURCE_MANIFEST.json。</p>'];
order=[numel(rec),1:numel(rec)-1];
for k=order,r=rec(k);stem=['figures/' r.name];html=[html '<section><h2>' r.family ' · ' r.name '</h2><p>' r.caption '</p><p><a href="' stem '.png">PNG</a><a href="' stem '.pdf">PDF</a><a href="' stem '.fig">MATLAB FIG</a></p><img loading="lazy" src="' r.png '"></section>'];end %#ok<AGROW>
html=[html '</html>'];writeutf(fullfile(root,'index.html'),html);
end
