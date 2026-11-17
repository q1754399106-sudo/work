function rerun_g1(names)
% Addendum (PLAN_STAGE2_ADD1_ZH.md): W1 and UQW1 for models whose stage-2 run preceded the addendum.
% Same data, noise estimate and t* fields as the main run (itag ''); output images2_<name>_g1.mat.
P=rtm_setup();
for n=names
 C=s2_cfg(n{1},1);C.tag='_g1';C.itag='';C.methods={'W1','UQW1'};C.fixed_methods=true;
 if ~isfile(fullfile(P.cache,['images2_' n{1} '_g1.mat'])),migrate2(C);end
 score_stage2(n{1},1);
end
end
