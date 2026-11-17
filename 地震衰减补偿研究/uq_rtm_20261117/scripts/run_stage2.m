function run_stage2(names,vscales)
% Stage-2 chain (PLAN_STAGE2_ZH.md) for the models in names (cell) and migration-velocity factors
% vscales (default 1): shots (if missing) -> datasets -> noise estimate -> migration -> scores.
% Every step is resumable / skipped when its output exists.
if nargin<2,vscales=1;end
P=rtm_setup();
for n=names
 name=n{1};fd_shots(name,1:numel(P.ixs));
 if ~isfile(fullfile(P.cache,['data_' name '.mat'])),prep_data(name);end
 for v=vscales
  C=s2_cfg(name,v);
  if isfile(fullfile(P.cache,['images2_' name C.tag '.mat'])),continue;end
  est_noise2(C);migrate2(C);score_stage2(name,v);
  if ~isempty(getenv('UQRTM_SMOKE')),return;end
 end
end
end
