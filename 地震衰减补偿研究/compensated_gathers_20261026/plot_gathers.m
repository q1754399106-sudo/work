function plot_gathers(force_models)
% All scientific graphics are rendered in MATLAB R2024b.
root=fileparts(mfilename('fullpath'));out=fullfile(root,'figures');if ~isfolder(out),mkdir(out);end
if nargin<1,force_models={};end
set(groot,'defaultFigureVisible','off','defaultAxesFontName','Microsoft YaHei','defaultTextFontName','Microsoft YaHei', ...
 'defaultAxesFontSize',10,'defaultTextInterpreter','none','defaultAxesTickLabelInterpreter','none','defaultLegendInterpreter','none');
diary(fullfile(root,'plot_MATLAB.log'));cl=onCleanup(@()diary('off'));maxNumCompThreads(4);
files=dir(fullfile(root,'results','*.mat'));records=struct('stem',{},'model',{},'method',{},'kind',{},'input_error',{},'output_error',{});
for k=1:numel(files)
 [~,stem]=fileparts(files(k).name);meta=fullfile(root,'results',[stem '.json']);if ~isfile(meta),continue;end
 r=jsondecode(fileread(meta));S=load(fullfile(files(k).folder,files(k).name));D=S.D;T=S.target;Y=S.compensated;t=S.t(:);E=Y-T;
 assert(isequal(size(D),size(T),size(Y)));assert(all(isfinite([D(:);T(:);Y(:)])));
 [x,xlab]=coords(r,size(D,2));mlabel=strrep(r.method,'+','＋');
 records(end+1)=rec(stem,r,'主图'); %#ok<AGROW>
 if ~triplet_exists(out,stem)||any(strcmp(r.model,force_models))
  lim=max(abs([D(:);T(:);Y(:);E(:)]));
  f=figure('Color','w','Units','pixels','Position',[10 10 1450 1190]);tl=tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact');
  ax=nexttile(tl);seismic(ax,x,t,T,lim,xlab);title(ax,'(a) 无衰减目标道集');
  ax=nexttile(tl);seismic(ax,x,t,D,lim,xlab);title(ax,sprintf('(b) 含衰减观测 · 14 dB · 误差 %.3f',r.input_error));
  ax=nexttile(tl);seismic(ax,x,t,Y,lim,xlab);title(ax,sprintf('(c) %s 补偿后 · 误差 %.3f',mlabel,r.output_error));
  ax=nexttile(tl);seismic(ax,x,t,E,lim,xlab);title(ax,'(d) 补偿后减目标（共用幅度色标）');
  live=find(any(D~=0,1));want=round([.33 .67]*size(D,2));
  for a=1:2
   [~,ii]=min(abs(live-want(a)));j=live(ii);ax=nexttile(tl);
   plot(ax,t,T(:,j),'Color',[.12 .15 .18],'LineWidth',1.25);hold(ax,'on');
   plot(ax,t,Y(:,j),'--','Color',[.8 .25 .08],'LineWidth',1.0);grid(ax,'on');xlim(ax,[t(1) t(end)]);
   xlabel(ax,'双程时间 (s)');ylabel(ax,'相对振幅');title(ax,sprintf('(%s) 第 %d 道 · 原始幅度波形',char('e'+a-1),j));legend(ax,{'目标','补偿'},'Location','best');
  end
  header(f,tl,[strrep(r.model,'_',' ') ' · ' mlabel ' 补偿道集'],context(r));
  deliver(f,out,stem);
 end
 if isfield(S,'QScan_Standard')||isfield(S,'Standard')
  cs=['compare_' stem];records(end+1)=rec(cs,r,'方法对照'); %#ok<AGROW>
  if ~triplet_exists(out,cs)||any(strcmp(r.model,force_models))
   A={T,D};labels={'无衰减目标','含衰减观测'};
   if isfield(S,'SC'),A{end+1}=S.SC;labels{end+1}='SC（未叠加 QScan）';end
   if isfield(S,'Standard'),A{end+1}=S.Standard;labels{end+1}='自研相干标准版';end
   if strcmp(r.method,'QScan+SC'),A{end+1}=Y;labels{end+1}='QScan＋SC';end
   if isfield(S,'QScan_Standard')&&~strcmp(r.model,'V4'),A{end+1}=S.QScan_Standard;labels{end+1}='QScan＋自研相干标准版';end
   lim=0;for j=1:numel(A),lim=max(lim,max(abs(A{j}),[],'all'));end
   rows=ceil(numel(A)/2);f=figure('Color','w','Units','pixels','Position',[10 10 1400 400*rows+100]);tl=tiledlayout(f,rows,2,'TileSpacing','compact','Padding','compact');
   for j=1:numel(A),ax=nexttile(tl);seismic(ax,x,t,A{j},lim,xlab);er=norm(A{j}-T,'fro')/norm(T,'fro');title(ax,sprintf('%s · 误差 %.3f',labels{j},er));end
   header(f,tl,[r.model ' · 原始候选方法道集对照'],{'全部面板共用幅度色标；方法参数固定，真值仅用于图中评分，不用于选择输出'});deliver(f,out,cs);
  end
 end
end
fid=fopen(fullfile(root,'FIGURE_INDEX.json'),'w','n','UTF-8');fwrite(fid,unicode2native(jsonencode(records,PrettyPrint=true),'UTF-8'));fclose(fid);
fprintf('PLOTTED %d figures from %d cases\n',numel(records),numel(files));
end

function r=rec(stem,m,kind)
r=struct('stem',stem,'model',m.model,'method',m.method,'kind',kind,'input_error',m.input_error,'output_error',m.output_error);
end
function [x,lab]=coords(r,n)
if ~isempty(regexp(r.model,'^H[1-6]$','once')),x=(0:n-1)*12;lab='横向距离 (m)';
elseif any(strcmp(r.model,{'V3','V5'})),x=(0:n-1)*30;lab='横向距离 (m)';
elseif strcmp(r.family,'四层模型'),x=(0:n-1)*50;lab='横向距离 (m)';
else,x=1:n;lab='道号';end
end
function seismic(ax,x,t,A,lim,xlab)
assert(max(abs(A),[],'all')<=lim+eps(lim)*8);
imagesc(ax,x,t,A,[-lim lim]);set(ax,'YDir','reverse','Layer','top');
colormap(ax,interp1([0 .5 1],[.10 .28 .50;1 1 1;.62 .13 .16],linspace(0,1,257)));
cb=colorbar(ax);cb.Label.String='相对振幅（统一色标）';xlabel(ax,xlab);ylabel(ax,'双程时间 (s)');
end
function s=context(r)
if contains(r.case,'Qminus30'),q=sprintf('Q 低估 30%%；QScan 因子 %.3f',r.qscan_factor);
elseif strcmp(r.method,'QScan+SC'),q=sprintf('四层标定 Q 不确定；QScan 因子 %.3f',r.qscan_factor);
else,q='给定参考 Q；SC 固定参数 a=0.125，μ相对强度=10';end
if r.solver.converged,c='已收敛';else,c='未达收敛容差';end
extra='';if ~isempty(r.dead_traces),extra=sprintf('；死道 %s（全道误差包含死道）',mat2str(r.dead_traces'));end
if strcmp(r.family,'四层模型'),extra=[extra '；新生成 FaultSection 展示观测'];end
if any(strcmp(r.model,{'F1','F2'})),extra=[extra '；SC 选参数据'];end
s={sprintf('%s；%s，KKT %.4f',q,c,r.solver.kkt_ratio),[extra '；无 AGC、无逐道归一化、无裁幅']};
end
function tf=triplet_exists(out,s)
tf=all(cellfun(@(ext)isfile(fullfile(out,[s ext])),{'.png','.pdf','.fig'}));
end
function header(f,tl,heading,sub)
% Reserve a fixed header strip: long subtitles must not push titles off canvas.
tl.Units='normalized';tl.OuterPosition=[0 .015 1 .87];
annotation(f,'textbox',[.015 .95 .97 .04],'String',heading,'LineStyle','none','FontSize',19,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none');
annotation(f,'textbox',[.015 .89 .97 .058],'String',sub,'LineStyle','none','FontSize',9.5,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none');
end
function deliver(f,out,stem)
drawnow;exportgraphics(f,fullfile(out,[stem '.png']),'Resolution',180,'BackgroundColor','white');
exportgraphics(f,fullfile(out,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');savefig(f,fullfile(out,[stem '.fig']),'compact');close(f);fprintf('Saved %s\n',stem);
end
