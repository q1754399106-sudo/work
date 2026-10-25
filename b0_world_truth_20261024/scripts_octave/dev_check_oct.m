function dev_check_oct(ids,tag)
% Stage 0 on Octave: diagnose2 on development gathers (already seen) with the same observations
% and bootstrap seeds as run_dev.m; output compared with the MATLAB analysis/dev_w*.csv.
P=oct_setup();M=read_csv_struct(fullfile(P.in,'qscan_bias_diagnostic_20261024','predictions','manifest.csv'));
out=fullfile(P.root,'octave_check',sprintf('dev_oct_%s.csv',tag));assert(~isfile(out),'exists');
rows=[];started=tic;
for id=ids
 g=M([M.id]==id);assert(g.q_err==0);[D,cm,Qused,offset]=observed_gather_oct(P,g);
 r=diagnose2(P,D,cm,Qused,offset,g.boot_seed);r.id=id;rows=[rows;r]; %#ok<AGROW>
 write_rows_csv(out,rows);
 fprintf('g%03d %s %s %s: c=%.4f B0 %+.6f B0M %+.6f B0C %+.6f (%.0fs, total %.0fs)\n',id,g.set,g.model,g.condition, ...
  r.qscan_factor,r.B0_rel_delta_mse,r.B0M_rel_delta_mse,r.B0C_rel_delta_mse,r.seconds,toc(started));
end
end
