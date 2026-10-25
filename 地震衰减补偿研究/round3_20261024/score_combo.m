function score_combo(cand)
% Truth-side scoring and PLAN_QSCAN_COH gates for the finished scenes.
% cand: combination under test (default 'QScan_Standard'; interim switch may use QScan_v1/QScan_Soft).
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));
if nargin<1,cand='QScan_Standard';end
methods={'FullCalRefit','SoftGroupRefit','Standard','QScanRefit','QScan_Soft','QScan_v1','QScan_Standard'};rows={};
for v={'V4','V5','F1','F2'}
 TT=load(fullfile(root,'dev',[v{1} '_truth.mat']));in=load(fullfile(root,'dev',[v{1} '.mat']),'config');
 for e=[-50 -30 -10 0 10 30 50]
  for k=1:2
   f=fullfile(root,'results',sprintf('%s_Q%+d_%s.mat',v{1},e,TT.tr(k).name));if ~isfile(f),continue;end
   r=load(f,'gathers','scan');T=struct('target',TT.tr(k).target,'masks',struct());
   if isfield(TT,'events'),T.reflectivity=TT.reflectivity;T.events=TT.events;T.weak=TT.weak;end
   if isfield(TT,'regions')
    if isfield(TT.regions,'fault'),T.masks.fault=TT.regions.fault;T.masks.pinch=TT.regions.pinch;end
    if isfield(TT.regions,'low_energy'),T.masks.low_energy=TT.regions.low_energy(:,:,k);end
   end
   for m=methods
    x=score_gather(r.gathers.(m{1}),T,in.config);
    rows{end+1}=struct('set',v{1},'condition',TT.tr(k).name,'Q_percent',e,'method',m{1},'factor',r.scan.factor, ...
     'overall',x.overall,'weak',x.weak,'fault',x.fault,'pinch',x.pinch,'low_energy',x.low_energy,'leak',x.leak); %#ok<AGROW>
   end
  end
 end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis','combo_metrics.csv'));
rat=@(sel,a,b,f)mean(M.(f)(sel&strcmp(M.method,a))./M.(f)(sel&strcmp(M.method,b)),'omitnan');
out=struct();
for grp={{'F1','F2'},{'V4','V5'},{'V4'},{'V5'},{'F1'},{'F2'}}
 sel=ismember(M.set,grp{1});key=strjoin(grp{1},'_');if ~any(sel),continue;end
 q0=sel&M.Q_percent==0;
 g=struct('vs_QScan_overall',rat(sel,cand,'QScanRefit','overall'),'vs_QScan_Q0',rat(q0,cand,'QScanRefit','overall'), ...
  'vs_QScan_weak',rat(sel,cand,'QScanRefit','weak'),'vs_QScan_fault',rat(sel,cand,'QScanRefit','fault'),'vs_QScan_pinch',rat(sel,cand,'QScanRefit','pinch'), ...
  'vs_Standard_Q0',rat(q0,cand,'Standard','overall'),'vs_Standard_worst',worst(M,sel,cand,'Standard'), ...
  'QScan_vs_FullCal_worst',worst(M,sel,'QScanRefit','FullCalRefit'),'scenes',sum(sel&strcmp(M.method,cand)));
 g.G1=g.vs_QScan_overall<=.97;g.G2=g.vs_QScan_Q0<=1.02;g.G3=g.vs_Standard_worst<=.95;g.G4=g.vs_Standard_Q0<=1.02;
 out.(key)=g;
 fprintf('%-6s n=%2d | vs QScan: all %.4f Q0 %.4f weak %.4f fault %.4f pinch %.4f | vs Standard: Q0 %.4f worst %.4f | QScan/FullCal worst %.4f | G1-4 %d%d%d%d\n', ...
  key,g.scenes,g.vs_QScan_overall,g.vs_QScan_Q0,g.vs_QScan_weak,g.vs_QScan_fault,g.vs_QScan_pinch,g.vs_Standard_Q0,g.vs_Standard_worst,g.QScan_vs_FullCal_worst,g.G1,g.G2,g.G3,g.G4);
end
uq_writejson(fullfile(root,'analysis',['combo_summary_' cand '.json']),struct('candidate',cand,'groups',out));
end

function w=worst(M,sel,a,b)
% Mean over (set, condition) of max-over-Q error ratio a/b.
keys=unique(strcat(M.set(sel),'_',M.condition(sel)));r=zeros(numel(keys),1);
for i=1:numel(keys)
 s=sel&strcmp(strcat(M.set,'_',M.condition),keys{i});r(i)=max(M.overall(s&strcmp(M.method,a)))/max(M.overall(s&strcmp(M.method,b)));
end
w=mean(r);
end
