function check_reproduction()
% Stage 0: the harness must regenerate the archived noise and sigma, reproduce the archived
% FullCalRefit/QScan outputs and old scores, and the split QScan (select + refit) must equal
% the original qscan_refit bit for bit. Also times one replicate.
P=bv_setup();maxNumCompThreads(P.threads);rows={};
ref=readtable(fullfile(P.repo,'theory_predict_diagnostic_20261024','analysis','case_scores.csv'),'TextType','string');
for s=1:numel(P.scenes)
 sc=load_scene(P,s);live=any(sc.D~=0,1);
 N=noise_replica(sc.D0,sc.snr_db,sc.archived_seed);
 noise_diff=max(abs(sc.D0+N-sc.D),[],'all')/max(abs(sc.D),[],'all');
 sigma_diff=max(abs(estimate_sigmas(sc.D,sc.cm.t)-sc.sigma))/max(sc.sigma);
 A=load(fullfile(P.bae,'results','ind',sprintf('%s_Q+0_%s.mat',sc.model,sc.condition)),'g');
 B=load(fullfile(P.bae,'results','ind_qscan',sprintf('%s_Q+0_%s.mat',sc.model,sc.condition)),'g','info');
 used=sc.cm;used.Q=sc.Q0;
 tic;[S1,i1]=full_refit(sc.D,sc.sigma,sc.cm,sc.Q0,sc.offset);t_refit=toc;
 tic;[Sq,iq]=qscan_refit(sc.D,sc.sigma,used,sc.offset);t_qscan=toc;
 tic;sel=qscan_select(sc.D,sc.sigma,sc.cm,sc.Q0,sc.offset);t_select=toc;
 Ss=full_refit(sc.D,sc.sigma,sc.cm,sc.Q0*P.grid(sel.best),sc.offset);
 base_diff=norm(S1-double(A.g.FullCalRefit),'fro')/norm(double(A.g.FullCalRefit),'fro');
 scan_diff=norm(Sq-double(B.g.QScanRefit),'fro')/norm(double(B.g.QScanRefit),'fro');
 split_diff=max(abs(Ss-Sq),[],'all');objective_diff=max(abs(sel.objective-iq.objective));
 terms_diff=max(abs(sel.misfit+sel.penalty-sel.objective)./abs(sel.objective));
 rel_base=norm(S1(:,live)-sc.T(:,live),'fro')/norm(sc.T(:,live),'fro');
 rel_scan=norm(Sq(:,live)-sc.T(:,live),'fro')/norm(sc.T(:,live),'fro');
 k=ref.model==sc.model&ref.condition==sc.condition;
 score_diff=max(abs([rel_base-ref.base_old_live(k),rel_scan-ref.scan_old_live(k)]));
 ok=noise_diff==0&&sigma_diff==0&&base_diff<1e-6&&scan_diff<1e-6&&split_diff==0&&objective_diff==0 ...
  &&iq.selected_factor==B.info.selected_factor&&terms_diff<1e-12&&score_diff<1e-6&&i1.optimization_flags==0;
 rows(end+1,:)={sc.model,sc.condition,noise_diff,sigma_diff,base_diff,scan_diff,split_diff,objective_diff,terms_diff, ...
  rel_base,rel_scan,score_diff,iq.selected_factor,B.info.selected_factor,i1.optimization_flags,t_refit,t_qscan,t_select,ok}; %#ok<AGROW>
 fprintf('%s %s noise %.3g sigma %.3g base %.3g scan %.3g split %.3g obj %.3g terms %.3g score %.3g factor %.4f ok=%d | refit %.1fs qscan %.1fs select %.1fs\n', ...
  sc.model,sc.condition,noise_diff,sigma_diff,base_diff,scan_diff,split_diff,objective_diff,terms_diff,score_diff,iq.selected_factor,ok,t_refit,t_qscan,t_select);
end
T=cell2table(rows,'VariableNames',{'model','condition','noise_replay_max_rel','sigma_replay_max_rel','base_vs_archived_rel', ...
 'qscan_vs_archived_rel','split_qscan_max_abs','objective_max_abs','terms_sum_rel','rel_base','rel_scan','old_score_max_abs', ...
 'selected_factor','archived_factor','optimization_flags','seconds_refit','seconds_qscan','seconds_select','passed'});
writetable(T,fullfile(P.root,'analysis','reproduction_check.csv'));
assert(all(T.passed),'Stage 0 reproduction failed.');
fprintf('DONE stage 0\n');
end
