function run_dev(worker,workers)
% Development worker: the 44 Q0 gathers of qscan_bias_diagnostic_20261024 (all already seen),
% same observations and bootstrap seeds. Writes truth-free values only; the original B0 must
% reproduce that round's sealed value.
P=bw_setup();maxNumCompThreads(P.threads);
M=readtable(fullfile(P.in,'qscan_bias_diagnostic_20261024','predictions','manifest.csv'),'TextType','string');
old=table();for w=1:4,old=[old;readtable(fullfile(P.in,'qscan_bias_diagnostic_20261024','predictions',sprintf('diagnostics_w%d.csv',w)))];end %#ok<AGROW>
M=M(M.q_err==0,:);M=M(mod((1:height(M))'-1,workers)==worker-1,:);
csv=fullfile(P.root,'analysis',sprintf('dev_w%d.csv',worker));assert(~isfile(csv),'Worker output exists; refusing to overwrite.');
rows=[];started=tic;
for i=1:height(M)
 g=M(i,:);[D,cm,Qused,offset]=observed_gather(P,g);r=diagnose2(P,D,cm,Qused,offset,g.boot_seed);
 r.id=g.id;r.sealed_B0=old.B0_rel_delta_mse(old.id==g.id);rows=[rows;struct2table(r)]; %#ok<AGROW>
 writetable(rows,csv);
 fprintf('g%03d %s %s %s: c=%.4f B0 %+.3f (sealed %+.3f) B0M %+.3f B0C %+.3f |m*|/|x1| %.3f |XC|/|x1| %.3f (%.0fs, total %.0fs)\n', ...
  g.id,g.set,g.model,g.condition,r.qscan_factor,r.B0_rel_delta_mse,r.sealed_B0,r.B0M_rel_delta_mse,r.B0C_rel_delta_mse, ...
  r.energy_mean_over_x1,r.energy_coh_over_x1,r.seconds,toc(started));
end
fprintf('DONE dev worker %d: %d gathers %.0fs\n',worker,height(M),toc(started));
end
