function va_run(names)
% Full chain of PLAN_ZH.md for the models in names (char or cell): migration (resumable), scores;
% the judgement once both models are scored. Console output is copied to logs/run_<names>.log.
if ischar(names),names={names};end
here=fileparts(mfilename('fullpath'));root=fileparts(here);
diary(fullfile(root,'logs',['run_' strjoin(names,'_') '.log']));diary on
fprintf('=== va_run %s start %s\n',strjoin(names,' '),datestr(now));
for n=names
 [~,P]=va_cfg(n{1});
 if ~isfile(fullfile(P.cache,['images_va_' n{1} '.mat'])),va_mig(n{1});end
 va_score(n{1});fprintf('scored %s %s\n',n{1},datestr(now));
end
if isfile(fullfile(root,'analysis','va_step1_D2.json'))&&isfile(fullfile(root,'analysis','va_step1_D4.json')),va_judge();end
fprintf('=== va_run %s done %s\n',strjoin(names,' '),datestr(now));diary off
end
