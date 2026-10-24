function check_inputs()
% Stage 0: (1) masked_refit with every row observed equals full_refit; (2) noise-free gathers
% and targets of all six geologies reproduce the archived observations with the archived seeds;
% (3) the noise-model window and AR(1) estimate per scene (truth-free statistics only).
P=bd_setup();maxNumCompThreads(P.threads);rows={};
sc=load_obs(P,'bae','F1','Moderate');S=full_refit(sc.D,sc.sigma,sc.cm,sc.Q0,sc.offset);
[Dp,Sm]=masked_refit(sc.D,sc.sigma,sc.cm,sc.Q0,sc.offset,true(sc.cm.nt,1));
c=sc.cm;c.Q=sc.Q0;d=prune_dictionary(c,sc.offset);
mask_diff=max(abs(Sm-S),[],'all');pred_diff=max(abs(Dp-d.G*S),[],'all')/max(abs(Dp),[],'all');
fprintf('masked_refit full mask: output diff %.3g, prediction vs G*S %.3g\n',mask_diff,pred_diff);
assert(mask_diff<1e-12&&pred_diff<1e-10,'masked_refit does not reproduce full_refit.');
seeds=struct('F1',2026102451,'F2',2026102461,'F3',2026102471,'F4',2026102481);
for set=["bae","qc","r3"]
 models=["F3","F4"];if set=="bae",models=["F1","F2"];end
 for model=models
  [D0,T]=regen_clean(P,char(set),char(model));
  switch set
   case 'bae',tt=load(fullfile(P.in,'bae_q_20261024','truth','independent',[char(model) '_truth.mat']));
   case 'qc',tt=load(fullfile(P.in,'qscan_criterion_20261024','truth','confirm',[char(model) '_truth.mat']));
   case 'r3',tt=load(fullfile(P.in,'round3_20261024','dev',[char(model) '_truth.mat']));
  end
  tdiff=max(abs(T-tt.tr(1).target),[],'all')/max(abs(T),[],'all');
  for k=1:2
   cond=["Moderate","Strong"];snr=[14 8];o=load_obs(P,char(set),char(model),char(cond(k)));
   N=noise_replica(D0,snr(k),seeds.(char(model))+k);ddiff=max(abs(D0+N-o.D),[],'all')/max(abs(o.D),[],'all');
   [phi,sigb,tend]=noise_model(o.D,o.cm);phi_old=noise_ar1(o.D,o.cm.t,P.window);
   rows(end+1,:)={char(set),char(model),char(cond(k)),ddiff,tdiff,tend,phi,phi_old,mean(sigb)/mean(o.sigma)}; %#ok<AGROW>
   fprintf('%s %s %s: obs replay %.3g, target %.3g, window end %.3f s, phi %.3f (full window %.3f), sigma ratio %.3f\n',rows{end,:});
  end
 end
end
T=cell2table(rows,'VariableNames',{'set','model','condition','observation_replay_rel','target_replay_rel','noise_window_end_s', ...
 'phi_hat','phi_hat_full_window','sigma_ratio_new_window'});
writetable(T,fullfile(P.root,'analysis','input_checks.csv'));
assert(all(T.observation_replay_rel<1e-12)&&all(T.target_replay_rel<1e-12),'Replay failed.');
fprintf('DONE stage 0\n');
end
