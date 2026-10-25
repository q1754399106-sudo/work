function score_sealed(root)
% Truth-side scoring of sealed validation and the frozen gates (Standard vs v1).
if nargin<1,root=fileparts(mfilename('fullpath'));end
base=fileparts(root);addpath(fullfile(root,'core'),fullfile(base,'coherence_soft_20261023','core'),fullfile(base,'coherence_soft_20261023','source'));
F=jsondecode(fileread(fullfile(root,'inputs','FROZEN.json')));rows={};
methods={'FullCalRefit','PrunedFullCalRefit','SoftGroupRefit','v1','Standard','Robust'};
for v={'V4','V5'}
 TT=load(fullfile(root,'truth',[v{1} '_truth.mat']));in=load(fullfile(root,'inputs','sealed',[v{1} '.mat']),'config');
 for e=[-30 -10 0 10 30]
  for k=1:2
   r=load(fullfile(root,'results',sprintf('%s_Q%+d_%s.mat',v{1},e,TT.tr(k).name)),'gathers');
   T=struct('target',TT.tr(k).target,'masks',struct('fault',TT.regions.fault,'pinch',TT.regions.pinch));
   if isfield(TT,'events'),T.reflectivity=TT.reflectivity;T.events=TT.events;T.weak=TT.weak;end
   if isfield(TT.regions,'low_energy'),T.masks.low_energy=TT.regions.low_energy(:,:,k);end
   for m=1:numel(methods)
    x=score_gather(r.gathers.(methods{m}),T,in.config);
    rows{end+1}=struct('set',v{1},'condition',TT.tr(k).name,'Q_percent',e,'method',methods{m},'overall',x.overall,'weak',x.weak, ...
     'fault',x.fault,'pinch',x.pinch,'low_energy',x.low_energy,'leak',x.leak); %#ok<AGROW>
   end
  end
 end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis','sealed_metrics.csv'));
ref=M(strcmp(M.method,'v1'),:);res=struct();
for m=methods
 x=M(strcmp(M.method,m{1}),:);rr=@(f,sel)mean(x.(f)(sel)./ref.(f)(sel),'omitnan');every=true(height(x),1);
 res.(m{1})=struct('overall',rr('overall',every),'weak',rr('weak',every),'fault',rr('fault',every),'pinch',rr('pinch',every), ...
  'low_energy',rr('low_energy',every),'leak',rr('leak',every),'Q_under',rr('overall',x.Q_percent<0), ...
  'V4_overall',rr('overall',strcmp(x.set,'V4')),'V5_overall',rr('overall',strcmp(x.set,'V5')), ...
  'correct_or_high_Q',rr('overall',x.Q_percent>=0),'scenes_gt_1',sum(x.overall>1),'mean_abs',mean(x.overall));
end
gt=F.sealed_validation.gates_for_Standard_vs_v1;s=res.Standard;
gates=struct('overall',s.overall<=gt.overall,'weak',s.weak<=gt.weak,'fault',s.fault<=gt.fault,'pinch',s.pinch<=gt.pinch,'Q_under',s.Q_under<=gt.Q_under_overall);
gf=fieldnames(gates);allpass=all(cellfun(@(f)logical(gates.(f)),gf)); % fixed after first run: a local variable named 'all' shadowed all()
uq_writejson(fullfile(root,'analysis','sealed_summary.json'),struct('ratios_vs_v1',res,'gates_Standard',gates,'all_pass',allpass));
fprintf('%-20s %7s %7s %7s %7s %7s %7s %7s %7s %7s %4s\n','method','overall','weak','fault','pinch','lowE','Qunder','Q>=0','V4','V5','>1');
for m=methods,r=res.(m{1});fprintf('%-20s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %4d\n',m{1},r.overall,r.weak,r.fault,r.pinch,r.low_energy,r.Q_under,r.correct_or_high_Q,r.V4_overall,r.V5_overall,r.scenes_gt_1);end
fprintf('Standard gates: %s all_pass=%d\n',jsonencode(gates),allpass);
end
