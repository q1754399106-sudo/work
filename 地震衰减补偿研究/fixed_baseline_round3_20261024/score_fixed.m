function score_fixed(sets)
% Score FixedSplit with round3's score_gather and the same truth/masks as
% score_guard.m. Also rescores round3's FullCalRefit on every scene and
% compares it with guard_metrics_*.csv as a reproduction check.
root=fileparts(mfilename('fullpath'));r3=fullfile(fileparts(root),'round3_20261024');
addpath(fullfile(r3,'core'),fullfile(r3,'source'));
if nargin<1,sets={'F3','F4','V5'};end
rows={};maxdiff=0;
for v=sets
 TT=load(fullfile(r3,'dev',[v{1} '_truth.mat']));in=load(fullfile(r3,'dev',[v{1} '.mat']),'config');
 ref=readtable(fullfile(r3,'analysis',sprintf('guard_metrics_%s.csv',guardname(v{1}))),'TextType','string');
 for e=[-50 -30 -10 0 10 30 50]
  for k=1:numel(TT.tr)
   f=fullfile(root,'results',sprintf('%s_Q%+d_%s.mat',v{1},e,TT.tr(k).name));
   T=struct('target',TT.tr(k).target,'masks',struct());
   if isfield(TT,'events'),T.reflectivity=TT.reflectivity;T.events=TT.events;T.weak=TT.weak;end
   if isfield(TT,'regions')&&isfield(TT.regions,'fault'),T.masks.fault=TT.regions.fault;T.masks.pinch=TT.regions.pinch;end
   r=load(f,'S');x=score_gather(r.S,T,in.config);
   rows{end+1}=struct('set',v{1},'condition',TT.tr(k).name,'Q_percent',e,'method','FixedSplit', ...
    'overall',x.overall,'weak',x.weak,'fault',x.fault,'pinch',x.pinch,'leak',x.leak); %#ok<AGROW>
   g=load(fullfile(r3,'results',sprintf('%s_Q%+d_%s.mat',v{1},e,TT.tr(k).name)),'gathers');
   y=score_gather(g.gathers.FullCalRefit,T,in.config);
   sel=ref.set==v{1}&ref.condition==TT.tr(k).name&ref.Q_percent==e&ref.method=="FullCalRefit";assert(nnz(sel)==1);
   maxdiff=max(maxdiff,abs(y.overall-ref.overall(sel))/ref.overall(sel));
  end
 end
end
M=struct2table([rows{:}]);writetable(M,fullfile(root,'analysis',sprintf('fixed_metrics_%s.csv',strjoin(sets,'_'))));
uq_writejson(fullfile(root,'analysis',sprintf('score_checks_%s.json',strjoin(sets,'_'))), ...
 struct('fullcal_rescore_max_relative_difference',maxdiff,'rows',height(M)));
fprintf('FullCalRefit rescore max rel diff %.3g; %d rows\n',maxdiff,height(M));
end

function n=guardname(s)
if any(strcmp(s,{'F3','F4'})),n='F3_F4';else,n=s;end
end
