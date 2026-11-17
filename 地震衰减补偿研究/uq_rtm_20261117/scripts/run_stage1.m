function run_stage1(name)
% Stage 1 chain for one model (no scores): shots -> eikonal fields -> datasets -> route C -> migration.
P=rtm_setup();maxNumCompThreads(P.threads);
fd_shots(name,1:numel(P.ixs));M=make_model(name,P);tstar_fields(M,P);
if ~isfile(fullfile(P.cache,['data_' name '.mat'])),prep_data(name);end
routec(name);est_noise(name);migrate(name);
end
