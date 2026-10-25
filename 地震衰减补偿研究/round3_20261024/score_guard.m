function score_guard(sets,tau)
% GuardedQScan (PLAN_ADDENDUM_GUARD.json): QScan route only when the relative
% objective gain exceeds tau, else the nominal route. Truth used for scoring only.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'core'),fullfile(root,'source'));
if nargin<1,sets={'F3','F4'};end;if nargin<2,tau=.01;end
rows={};
for v=sets
 TT=load(fullfile(root,'dev',[v{1} '_truth.mat']));in=load(fullfile(root,'dev',[v{1} '.mat']),'config');
 for e=[-50 -30 -10 0 10 30 50]
  for k=1:2
   f=fullfile(root,'results',sprintf('%s_Q%+d_%s.mat',v{1},e,TT.tr(k).name));if ~isfile(f),continue;end
   r=load(f,'gathers','scan');[~,k1]=min(abs(log(r.scan.grid)));gain=(r.scan.objective(k1)-min(r.scan.objective))/r.scan.objective(k1);use=gain>tau;
   g=r.gathers;
   if use,g.Guarded_Soft=g.QScan_Soft;g.Guarded_Standard=g.QScan_Standard;g.GuardedQScanRefit=g.QScanRefit;
   else,g.Guarded_Soft=g.SoftGroupRefit;g.Guarded_Standard=g.Standard;g.GuardedQScanRefit=g.FullCalRefit;end
   T=struct('target',TT.tr(k).target,'masks',struct());
   if isfield(TT,'events'),T.reflectivity=TT.reflectivity;T.events=TT.events;T.weak=TT.weak;end
   if isfield(TT,'regions')&&isfield(TT.regions,'fault'),T.masks.fault=TT.regions.fault;T.masks.pinch=TT.regions.pinch;end
   for m=fieldnames(g)'
    x=score_gather(g.(m{1}),T,in.config);
    rows{end+1}=struct('set',v{1},'condition',TT.tr(k).name,'Q_percent',e,'method',m{1},'qscan_used',use,'gain',gain, ...
     'factor',r.scan.factor,'overall',x.overall,'weak',x.weak,'fault',x.fault,'pinch',x.pinch); %#ok<AGROW>
   end
  end
 end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis',sprintf('guard_metrics_%s.csv',strjoin(sets,'_'))));
rat=@(sel,a,b)mean(M.overall(sel&strcmp(M.method,a))./M.overall(sel&strcmp(M.method,b)));
sel=true(height(M),1);q0=M.Q_percent==0;res=struct();
for c={{'Guarded_Soft','SoftGroupRefit'},{'Guarded_Standard','Standard'},{'QScan_Soft','SoftGroupRefit'},{'QScan_Standard','Standard'},{'GuardedQScanRefit','FullCalRefit'},{'QScanRefit','FullCalRefit'}}
 a=c{1}{1};b=c{1}{2};
 x=struct('H1_vs_QScanRefit',rat(sel,a,'QScanRefit'),'H2_Q0_vs_FullCal',rat(q0,a,'FullCalRefit'), ...
  'H3_worst_vs_FullCal',worst(M,a,'FullCalRefit'),'H4_Q0_vs_own_nominal',rat(q0,a,b),'mean_overall',mean(M.overall(strcmp(M.method,a))), ...
  'weak_vs_QScanRefit',mean(M.weak(strcmp(M.method,a))./M.weak(strcmp(M.method,'QScanRefit')),'omitnan'));
 x.pass=x.H1_vs_QScanRefit<=.97&&x.H2_Q0_vs_FullCal<=1&&x.H3_worst_vs_FullCal<=.5&&x.H4_Q0_vs_own_nominal<=1.02;res.(a)=x;
 fprintf('%-18s H1 %.4f  H2 %.4f  H3 %.4f  H4 %.4f  mean %.4f  weak %.4f  pass %d\n',a,x.H1_vs_QScanRefit,x.H2_Q0_vs_FullCal,x.H3_worst_vs_FullCal,x.H4_Q0_vs_own_nominal,x.mean_overall,x.weak_vs_QScanRefit,x.pass);
end
used=M(strcmp(M.method,'QScanRefit'),:);fprintf('QScan applied in %d/%d scenes (Q0: %d/%d)\n',sum(used.qscan_used),height(used),sum(used.qscan_used&used.Q_percent==0),sum(used.Q_percent==0));
uq_writejson(fullfile(root,'analysis',sprintf('guard_summary_%s.json',strjoin(sets,'_'))),struct('tau',tau,'sets',{sets},'results',res));
end

function w=worst(M,a,b)
keys=unique(strcat(M.set,'_',M.condition));r=zeros(numel(keys),1);
for i=1:numel(keys)
 s=strcmp(strcat(M.set,'_',M.condition),keys{i});r(i)=max(M.overall(s&strcmp(M.method,a)))/max(M.overall(s&strcmp(M.method,b)));
end
w=mean(r);
end
